import os
import pandas as pd

def load_data(files):
    data_dict = {}
    for file in files:
        name = file.split('.')[0]
        data_dict[name] = pd.read_csv(file, sep='\t', index_col=0)
    return data_dict

def get_genes(data_dict):
    genes_dict = {}
    for cancer, data in data_dict.items():
        genes_dict[cancer] = set(data.columns)
    return genes_dict

def get_common_genes(genes_dict):
    return set.intersection(*genes_dict.values())

def get_specific_or_overlap(cancer, genes, genes_dict, common_genes):
    other_genes = set.union(*(genes_dict[other] for other in genes_dict if other != cancer))
    specific_genes = genes - other_genes
    if specific_genes:
        return 'specific', specific_genes
    else:
        overlap_genes = genes - common_genes
        return 'overlapping', overlap_genes

def save_gene_results(common_genes, specific_or_overlap_genes, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    with open(os.path.join(output_dir, 'common_genes.txt'), 'w') as f:
        f.write('\n'.join(common_genes))

    for cancer, (gene_type, genes) in specific_or_overlap_genes.items():
        file_path = os.path.join(output_dir, f'{cancer}_{gene_type}_genes.txt')
        with open(file_path, 'w') as f:
            f.write('\n'.join(genes))

def main():
    cancer_files = [
        'breast cancer (BRACA).tsv', 'Adenoid cystic carcinoma (ACC).tsv', 'Bladder urothelial carcinoma (BLCA).tsv',
        'Cervical Squamous Cell Carcinoma and Endocervical Adenocarcinoma (CESE).tsv', 'Cholangiocarcinoma (CHOL).tsv',
        'Diffuse large B-cell lymphoma (DLBCL).tsv', 'Esophageal cancer (ESCA).tsv', 'glioblastoma multiforme (GBM).tsv',
        'Head and neck squamous cell carcinoma (HNSC).tsv', 'Kidney renal clear cell carcinoma (KIRC).tsv',
        'Kidney renal papillary cell carcinoma (KIRP).tsv', 'Acute myeloid leukemia (LAML).tsv',
        'Pediatric low-grade gliomas (LGG).tsv', 'Liver Hepatocellular Carcinoma(LIHC).tsv', 'Lung adenocarcinoma (LUAD).tsv',
        'lung squamous cell carcinoma (LUSC).tsv', 'mesothelioma cancer (MESO).tsv', 'Ovarian cancer (OV).tsv',
        'Pancreatic adenocarcinoma (PAAD).tsv', 'Pheochromocytoma and paraganglioma (PCPGs).tsv',
        'Prostate Adenocarcinoma (PRAD).tsv', 'Rectum Adenocarcinoma (READ).tsv', 'Sarcoma Cancer (SARC).tsv',
        'Skin cutaneous melanoma (SKCM).tsv', 'Stomach adenocarcinoma (STAD).tsv', 'Tenosynovial giant cell tumors (TGCTs).tsv',
        'Thyroid carcinoma (THCA).tsv', 'Thymoma and Thymic Carcinoma(THYM).tsv', 'Uterine corpus endometrial carcinoma (UCEC).tsv',
        'Uterine Carcinosarcoma (UCS).tsv', 'Uveal melanoma (UVM).tsv', 'Colon adenocarcinoma (COAD).tsv',
        'Kidney Chromophobe (KICH).tsv'
    ]

    cancer_data = load_data(cancer_files)
    gene_sets = get_genes(cancer_data)
    common_genes = get_common_genes(gene_sets)
    print("Common Genes:", common_genes)

    specific_or_overlap_genes = {}
    for cancer, genes in gene_sets.items():
        gene_type, found_genes = get_specific_or_overlap(cancer, genes, gene_sets, common_genes)
        specific_or_overlap_genes[cancer] = (gene_type, found_genes)
        print(f"\nGenes for {cancer}: {len(found_genes)}")
        print(found_genes)

    output_folder = 'gene_results'
    save_gene_results(common_genes, specific_or_overlap_genes, output_folder)

if __name__=="__main__":
    main()