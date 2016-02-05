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
