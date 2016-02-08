from gemini import GeminiQuery


def var_is_rare(variant_data, threshold):
    """Check if variant is rare, as defined by the passed cutoff
    :param variant_data: A GeminiRow for a single variant.
    :type variant_data: GeminiRow.
    :param threshold: Allele frequency rarity threshold.
    :type threshold: float.
    :returns:  bool -- True or False.
    """

    if variant_data['in_esp'] != 0 or variant_data['in_1kg'] != 0 or variant_data['in_exac'] != 0:
        if variant_data['max_aaf_all'] > threshold:
            return False
        else:
            return True
    else:
        return True


def var_is_in_cosmic(variant_data):
    """Check if variant is in the COSMIC database
    :param variant_data: A GeminiRow for a single variant.
    :type variant_data: GeminiRow.
    :returns:  bool -- True or False.
    """

    if variant_data['cosmic_ids'] is not None:
        return True
    else:
        return False


def var_is_in_clinvar(variant_data):
    """Check if variant is in the ClinVar database
    :param variant_data: A GeminiRow for a single variant.
    :type variant_data: GeminiRow.
    :returns:  bool -- True or False.
    """

    if variant_data['clinvar_sig'] is not None:
        return True
    else:
        return False


def var_is_pathogenic(variant_data):
    """Check if variant is listed as pathogenic in ClinVar
    :param variant_data: A GeminiRow for a single variant.
    :type variant_data: GeminiRow.
    :returns:  bool -- True or False.
    """

    if variant_data['clinvar_sig'] is not None:
        if "pathogenic" in variant_data['clinvar_sig']:
            return True
        else:
            return False
    else:
        return False


def var_is_protein_effecting(variant_data):
    """Check if variant has a MED or HIGH impact
    :param variant_data: A GeminiRow for a single variant.
    :type variant_data: GeminiRow.
    :returns:  bool -- True or False.
    """
    if variant_data['impact_severity'] != "LOW":
        return True
    else:
        return False


def var_in_gene(variant_data, genes):
    """Check if variant has a gene name associated with it
    :param variant_data: A GeminiRow for a single variant.
    :type variant_data: GeminiRow.
    :returns:  bool -- True or False.
    """
    if variant_data['gene'] in genes:
        return True
    else:
        return False


def var_is_lof(variant_data):
    if variant_data['is_lof']:
        return True
    else:
        return False


def var_is_coding(variant_data):
    if variant_data['is_coding']:
        return True
    else:
        return False


def var_is_splicing(variant_data):
    if variant_data['is_splicing']:
        return True
    else:
        return False


def parse_rs_ids(variant_data):
    if variant_data['rs_ids'] is not None:
        return variant_data['rs_ids'].split(',')
    else:
        return []


def parse_cosmic_ids(variant_data):
    if variant_data['cosmic_ids'] is not None:
        return variant_data['cosmic_ids'].split(',')
    else:
        return []


def gemini_query(db):
    """Execute a GEMINI Query
    :param db: A GEMINI database name
    :type db: str
    :returns:  GeminiRow Results -- True or False.
    """
    query = "SELECT chrom, start, end, ref, alt, vcf_id, rs_ids, cosmic_ids, filter, qual, qual_depth, depth, " \
            "type, sub_type, " \
            "gene, transcript, exon, codon_change, aa_change, biotype, impact, impact_so, impact_severity, " \
            "aa_length, is_lof, is_conserved, pfam_domain, in_omim, clinvar_sig, clinvar_disease_name, " \
            "is_exonic, is_coding, is_splicing, " \
            "clinvar_origin, clinvar_causal_allele, clinvar_dbsource, clinvar_dbsource_id, " \
            "clinvar_on_diag_assay, rmsk, in_segdup, strand_bias, rms_map_qual, in_hom_run, num_mapq_zero, " \
            "num_reads_w_dels, grc, gms_illumina, in_cse, num_alleles, allele_count, haplotype_score, " \
            "is_somatic, somatic_score, aaf_esp_ea, aaf_esp_aa, aaf_esp_all, aaf_1kg_amr, " \
            "aaf_1kg_eas, aaf_1kg_sas, aaf_1kg_afr, aaf_1kg_eur, aaf_1kg_all, aaf_exac_all, aaf_adj_exac_all, " \
            "aaf_adj_exac_afr, aaf_adj_exac_amr, aaf_adj_exac_eas, aaf_adj_exac_fin, aaf_adj_exac_nfe, " \
            "aaf_adj_exac_oth, aaf_adj_exac_sas, max_aaf_all, in_esp, in_1kg, in_exac, info," \
            "(gts).(*), (gt_depths).(*), (gt_ref_depths).(*), (gt_alt_depths).(*) FROM variants"

    gq = GeminiQuery(db)
    gq.run(query)

    return gq
