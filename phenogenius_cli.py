import numpy as np
import pandas as pd
import ujson as json
import pickle as pk
import sklearn
from pandarallel import pandarallel
import sys
import click
import logging


def load_data():
    matrix = pd.read_csv(
        "data/resources/ohe_all_thesaurus_weighted.tsv.gz",
        sep="\t",
        compression="gzip",
        index_col=0,
        engine="pyarrow",
    )
    return matrix


def load_nmf_model():
    with open("data/resources/pheno_NMF_390_model_42.pkl", "rb") as pickle_file:
        pheno_NMF = pk.load(pickle_file)
    with open("data/resources/pheno_NMF_390_matrix_42.pkl", "rb") as pickle_file:
        reduced = pk.load(pickle_file)
    return pheno_NMF, reduced


def symbol_to_id_to_dict():
    # from NCBI
    ncbi_df = pd.read_csv("data/resources/Homo_sapiens.gene_info.gz", sep="\t")
    ncbi_df = ncbi_df[ncbi_df["#tax_id"] == 9606]
    ncbi_df_ncbi = ncbi_df.set_index("Symbol")
    ncbi_to_dict_ncbi = ncbi_df_ncbi["GeneID"].to_dict()
    ncbi_df = ncbi_df.set_index("GeneID")
    ncbi_to_dict = ncbi_df["Symbol"].to_dict()
    return ncbi_to_dict_ncbi, ncbi_to_dict


def load_similarity_dict():
    with open("data/resources/similarity_dict_threshold_80.json") as json_data:
        data_dict = json.load(json_data)
    return data_dict


def load_hp_ontology():
    with open("data/resources/hpo_obo.json") as json_data:
        data_dict = json.load(json_data)
    return data_dict


def get_symbol(gene, symbol):
    if gene in symbol.keys():
        return symbol[gene]


def get_similar_terms(hpo_list, similarity_terms_dict):
    hpo_list_w_simi = {}
    for term in hpo_list:
        hpo_list_w_simi[term] = 1
        if term in similarity_terms_dict.keys():
            for key, value in similarity_terms_dict[term].items():
                if value > 0.8:
                    score = value / len(similarity_terms_dict[term].keys())
                    if key in hpo_list_w_simi.keys():
                        if score > hpo_list_w_simi[key]:
                            hpo_list_w_simi[key] = score
                            # print(key)
                            # print(statistics.mean(similarity_terms_dict[term].values()))
                        else:
                            pass
                    else:
                        hpo_list_w_simi[key] = score
    hpo_list_all = hpo_list_w_simi.keys()
    return hpo_list_w_simi, list(hpo_list_all)


def score_sim_add(hpo_list_add, matrix, sim_dict, symbol):
    matrix_filter = matrix[hpo_list_add]
    for key, value in sim_dict.items():
        matrix_filter[key] = matrix_filter[key] * value
    matrix_filter["sum"] = matrix_filter.sum(axis=1)
    matrix_filter["gene_symbol"] = matrix_filter.index.to_series().apply(
        get_symbol, args=(symbol,)
    )
    return matrix_filter.sort_values("sum", ascending=False)


@click.command()
@click.option("--result_file", default="match.tsv")
@click.option("--hpo_list", default="HP:0000107,HP:0000108,HP:0001407,HP:0005562")
def evaluate_matching(result_file, hpo_list):
    logging.info("INFO: load databases")
    ncbi, symbol = symbol_to_id_to_dict()
    data = load_data()
    hp_onto = load_hp_ontology()

    logging.info("INFO: clean HPO list")
    hpo_list_ini = hpo_list.strip().split(",")
    hpo_list_up = []

    for hpo in hpo_list_ini:
        if hpo in ["HP:0000001"]:
            pass
        else:
            if data[hpo].astype(bool).sum(axis=0) != 0:
                hpo_list_up.append(hpo)
            else:
                hpo_to_test = hp_onto[hpo]["direct_parent"][0]
                while data[hpo_to_test].astype(bool).sum(
                    axis=0
                ) == 0 and hpo_to_test not in ["HP:0000001"]:
                    hpo_to_test = hp_onto[hpo_to_test]["direct_parent"][0]
                if hpo_to_test in ["HP:0000001"]:
                    pass
                else:
                    hpo_list_up.append(hpo_to_test)
    hpo_list = list(set(hpo_list_up))
    if len(hpo_list) < 6:
        logging.info("INFO: selected symptom interaction model - NMF")
        pandarallel.initialize(nb_workers=4)
        pheno_NMF, reduced = load_nmf_model()
        witness = np.zeros(len(data.columns))
        witness_nmf = np.matmul(pheno_NMF.components_, witness)
        witness_df = (
            pd.DataFrame(reduced)
            .set_index(data.index)
            .parallel_apply(lambda x: sum((x - witness_nmf) ** 2), axis=1)
        )

        patient = np.zeros(len(data.columns))
        for hpo in hpo_list:
            hpo_index = list(data.columns).index(hpo)
            patient[hpo_index] = 1

        patient_nmf = np.matmul(pheno_NMF.components_, patient)
        patient_df = (
            pd.DataFrame(reduced)
            .set_index(data.index)
            .parallel_apply(lambda x: sum((x - patient_nmf) ** 2), axis=1)
        )

        case_df = pd.DataFrame(patient_df - witness_df)
        case_df.columns = ["score"]
        case_df["score_norm"] = abs(case_df["score"] - case_df["score"].max())
        case_df["sum"] = case_df["score_norm"]
        case_df_sort = case_df.sort_values(by="sum", ascending=False)
        case_df_sort["rank"] = (
            case_df_sort["sum"].rank(ascending=False, method="max").astype(int)
        )
        case_df_sort["gene_symbol"] = case_df_sort.index.to_series().apply(
            get_symbol, args=(symbol,)
        )
        match_nmf = case_df_sort[["gene_symbol", "rank", "sum"]]
        match_nmf_filter = match_nmf[match_nmf["sum"] > 0.01].reset_index()
        match_nmf_filter.columns = ["gene_id", "gene_symbol", "rank", "sum"]
        match_nmf_filter["sum"] = match_nmf_filter["sum"].round(2)
        match_nmf_filter.to_csv(result_file, sep="\t", index=False)

    else:
        logging.info("INFO: selected symptom interaction model - node similarity")
        similarity_terms_dict = load_similarity_dict()
        sim_dict, hpo_list_add = get_similar_terms(hpo_list, similarity_terms_dict)
        results_sum_add = score_sim_add(hpo_list_add, data, sim_dict, symbol)
        results_sum_add["rank"] = (
            results_sum_add["sum"].rank(ascending=False, method="max").astype(int)
        )
        cols = results_sum_add.columns.tolist()
        cols = cols[-2:] + cols[:-2]
        match_sim = results_sum_add[cols].sort_values(by=["sum"], ascending=False)
        match_sim_filter = match_sim[match_sim["sum"] > 0.01].reset_index()
        match_sim_filter_print = match_sim_filter.iloc[:, [0, 1, 2, -1]]
        match_sim_filter_print.columns = ["gene_id", "gene_symbol", "rank", "sum"]
        match_sim_filter_print["sum"] = match_sim_filter_print["sum"].round(2)
        match_sim_filter_print.to_csv(result_file, sep="\t", index=False)


if __name__ == "__main__":
    evaluate_matching()
