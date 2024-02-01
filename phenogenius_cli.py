import numpy as np
import pandas as pd
import ujson as json
import pickle as pk
from pandarallel import pandarallel
import click
import logging
import subprocess

def get_version_from_git_tag():
    try:
        version = subprocess.check_output(["git", "describe", "--tags"]).strip().decode('utf-8')
    except Exception:
        version = "Unknown"
    return version

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

def get_phenotype_specificity(row):
    rank = row["rank"]
    score = row["score"]
    if score < 0.1:
        return "D - the reported phenotype is NOT consistent with what is expected for the gene/genomic region or not consistent in general."
    elif rank < 41:
        return "A - the reported phenotype is highly specific and relatively unique to the gene (top 40, 50 perc of diagnosis in PhenoGenius cohort)."
    elif rank < 250:
        return "B - the reported phenotype is consistent with the gene, is highly specific, but not necessarily unique to the gene (top 250, 75 perc of diagnosis in PhenoGenius cohort)."
    else:
        return "C - the phenotype is reported with limited association with the gene, not highly specific and/or with high genetic heterogeneity."
    
def get_hpo_implicated_dict(data, hpo_list, hp_onto):
    data_filter = data[hpo_list]
    data_filter_dict = data_filter.to_dict(orient='index')
    annot_dict = {}
    for key, value in data_filter_dict.items():
        hpo_implicated = []
        hpo_description_implicated = []
        for k, v in value.items():
            if v > 0:
                hpo_implicated.append({k:round(v,1)})
                hpo_description_implicated.append({hp_onto[k]['name']:round(v,1)})

        annot_dict[key] = {"hpo_implicated": hpo_implicated, "hpo_description_implicated": hpo_description_implicated}
    return annot_dict

def add_hpo_implicated(x, annot_dict):
    if x in annot_dict.keys():
        return annot_dict[x]["hpo_implicated"]
    else:
        return None

def add_hpo_description_implicated(x, annot_dict):
    if x in annot_dict.keys():
        return annot_dict[x]["hpo_description_implicated"]
    else:
        return None
    
@click.command()
@click.version_option(version=get_version_from_git_tag(), prog_name="PhenoGenius")
@click.option("--result_file", default="match.tsv", help="Output file name, default = match.tsv")
@click.option("--hpo_list", default=None, help="(Mandatory) List of HPO terms to match, separated with commas")
@click.option("--gene_list", default=None, help="(Optional) List of genes in NCBI ID format to match, separated with commas")
def evaluate_matching(result_file, hpo_list, gene_list):
    logging.info("INFO: load databases")
    ncbi, symbol = symbol_to_id_to_dict()
    data = load_data()
    hp_onto = load_hp_ontology()

    logging.info("INFO: clean HPO list")
    if hpo_list is None:
        print("Provide a list of HPO terms with the --hpo_list option")
    else:
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

        annot_dict = get_hpo_implicated_dict(data, hpo_list, hp_onto)

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
            match_nmf_filter.columns = ["#gene_id", "gene_symbol", "rank", "score"]
            match_nmf_filter["score"] = match_nmf_filter["score"].round(2)
            match_nmf_filter["hpo_implicated"] = match_nmf_filter["#gene_id"].apply(add_hpo_implicated, args=(annot_dict,))
            match_nmf_filter["hpo_description_implicated"] = match_nmf_filter["#gene_id"].apply(add_hpo_description_implicated, args=(annot_dict,))
            match_nmf_filter["phenotype_specificity"] = match_nmf_filter.apply(get_phenotype_specificity, axis=1)
            if gene_list is not None:
                gene_list = gene_list.strip().split(",")
                gene_list_int = [eval(x) for x in gene_list]
                match_nmf_filter = match_nmf_filter[match_nmf_filter["#gene_id"].isin(gene_list_int)]
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
            match_sim_filter_print.columns = ["#gene_id", "gene_symbol", "rank", "score"]
            match_sim_filter_print["score"] = match_sim_filter_print["score"].round(2)
            match_sim_filter_print["hpo_implicated"] = match_sim_filter_print["#gene_id"].apply(add_hpo_implicated, args=(annot_dict,))
            match_sim_filter_print["hpo_description_implicated"] = match_sim_filter_print["#gene_id"].apply(add_hpo_description_implicated, args=(annot_dict,))
            match_sim_filter_print["phenotype_specificity"] = match_sim_filter_print.apply(get_phenotype_specificity, axis=1)
            if gene_list is not None:
                gene_list = gene_list.strip().split(",")
                gene_list_int = [eval(x) for x in gene_list]
                match_sim_filter_print = match_sim_filter_print[match_sim_filter_print["#gene_id"].isin([eval(i) for i in gene_list_int])]
            match_sim_filter_print.to_csv(result_file, sep="\t", index=False)


if __name__ == "__main__":
    evaluate_matching()
