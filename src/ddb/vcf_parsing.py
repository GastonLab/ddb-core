import sys
from cyvcf2 import VCF
from collections import defaultdict


def parse_caller_vcfs(sample_dict, caller_list):
    caller_vcf_records = defaultdict(lambda: dict())
    for caller in caller_list:
        parse_vcf(sample_dict[caller], caller, caller_vcf_records)

    return caller_vcf_records


def parse_vcf(vcf_file, caller, caller_vcf_records):
    sys.stdout.write("Reading {}\n".format(vcf_file))
    vcf = VCF(vcf_file)
    for record in vcf:
        if len(record.ALT) > 1:
            sys.stderr.write("ERROR: More than one alternative allele detected in file "
                             "{}\n Record: {}\n".format(vcf_file, record))
            sys.exit()
        key = (unicode("chr{}".format(record.CHROM)), int(record.start), int(record.end), unicode(record.REF),
               unicode(record.ALT[0]))
        caller_vcf_records[caller][key] = record


def parse_mutect_vcf_record(record):
    # Pseudocount. Encountered a division by zero issue in at least one mutect record
    depth = int(record.gt_depths[0])
    if depth < 1:
        depth = 1
    info = {'DP': str(depth),
            'FILTER': str(record.FILTER),
            'GTF_DP': str(record.gt_depths[0]),
            'GTF_AD': str(record.gt_alt_depths[0]),
            'MULTIALLELIC': str(record.INFO.get('OLD_MULTIALLELIC')) or None,
            'AAF': str(float(record.gt_alt_depths[0]) / float(depth))}

    return info


def parse_vardict_vcf_record(record):
    info = {'DP': str(record.INFO.get('DP')),
            'VD': str(record.INFO.get('VD')),
            'AF': str(record.INFO.get('AF')),
            'FILTER': str(record.FILTER),
            'BIAS': str(record.INFO.get('BIAS')),
            'REFBIAS': str(record.INFO.get('REFBIAS')),
            'VARBIAS': str(record.INFO.get('VARBIAS')),
            'QUAL': str(record.INFO.get('QUAL')),
            'QSTD': str(record.INFO.get('QSTD')),
            'SBF': str(record.INFO.get('SBF')),
            'ODDRATIO': str(record.INFO.get('ODDRATIO')),
            'MQ': str(record.INFO.get('MQ')),
            'SN': str(record.INFO.get('SN')),
            'HIAF': str(record.INFO.get('HIAF')),
            'ADJAF': str(record.INFO.get('ADJAF')),
            'MSI': str(record.INFO.get('MSI')),
            'MSILEN': str(record.INFO.get('MSILEN')),
            'SHIFT3': str(record.INFO.get('SHIFT3')),
            'NM': str(record.INFO.get('NM')),
            'GDAMP': str(record.INFO.get('GDAMP')),
            'LSEQ': str(record.INFO.get('LSEQ')),
            'RSEQ': str(record.INFO.get('RSEQ')),
            'TLAMP': str(record.INFO.get('TLAMP')),
            'NCAMP': str(record.INFO.get('NCAMP')),
            'AMPFLAG': str(record.INFO.get('AMPFLAG')),
            'HICNT': str(record.INFO.get('HICNT')),
            'HICOV': str(record.INFO.get('HICOV')),
            'GTF_DP': str(record.gt_depths[0]),
            'GTF_AD': str(record.gt_alt_depths[0]),
            'MULTIALLELIC': str(record.INFO.get('OLD_MULTIALLELIC')) or None,
            'AAF': str(float(record.gt_alt_depths[0]) / float(record.gt_depths[0]))}

    return info


def parse_freebayes_vcf_record(record):
    info = {'DP': str(record.INFO.get('DP')),
            'AF': str(record.INFO.get('AF')),
            'FILTER': str(record.FILTER),
            'AC': str(record.INFO.get('AC')),
            'RO': str(record.INFO.get('RO')),
            'AO': str(record.INFO.get('AO')),
            'PRO': str(record.INFO.get('PRO')),
            'PAO': str(record.INFO.get('PAO')),
            'QR': str(record.INFO.get('QR')),
            'QA': str(record.INFO.get('QA')),
            'PQR': str(record.INFO.get('PQR')),
            'PQA': str(record.INFO.get('PQA')),
            'SRF': str(record.INFO.get('SRF')),
            'SRR': str(record.INFO.get('SRR')),
            'SAF': str(record.INFO.get('SAF')),
            'SAR': str(record.INFO.get('SAR')),
            'SRP': str(record.INFO.get('SRP')),
            'SAP': str(record.INFO.get('SAP')),
            'AB': str(record.INFO.get('AB')),
            'ABP': str(record.INFO.get('ABP')),
            'RUN': str(record.INFO.get('RUN')),
            'RPP': str(record.INFO.get('RPP')),
            'RPPR': str(record.INFO.get('RPPR')),
            'RPL': str(record.INFO.get('RPL')),
            'RPR': str(record.INFO.get('RPR')),
            'EPP': str(record.INFO.get('EPP')),
            'EPPR': str(record.INFO.get('EPPR')),
            'DRPA': str(record.INFO.get('DRPA')),
            'ODDS': str(record.INFO.get('ODDS')),
            'GTI': str(record.INFO.get('GTI')),
            'TYPE': str(record.INFO.get('TYPE')),
            'CIGAR': str(record.INFO.get('CIGAR')),
            'NUMALT': str(record.INFO.get('NUMALT')),
            'MEANALT': str(record.INFO.get('MEANALT')),
            'LEN': str(record.INFO.get('LEN')),
            'MQM': str(record.INFO.get('MQM')),
            'MQMR': str(record.INFO.get('MQMR')),
            'PAIRED': str(record.INFO.get('PAIRED')),
            'PAIREDR': str(record.INFO.get('PAIREDR')),
            'GTF_DP': str(record.gt_depths[0]),
            'MULTIALLELIC': str(record.INFO.get('OLD_MULTIALLELIC')) or None,
            'AAF': str(float(record.INFO.get('AO')) / float(record.gt_depths[0]))}

    return info


def parse_scalpel_vcf_record(record):
    info = {'DP': str(record.gt_depths[0]),
            'AVGCOV': str(record.INFO.get('AVGCOV')),
            'MINCOV': str(record.INFO.get('MINCOV')),
            'ALTCOV': str(record.INFO.get('ALTCOV')),
            'COVRATIO': str(record.INFO.get('COVRATIO')),
            'FILTER': str(record.FILTER),
            'ZYG': str(record.INFO.get('ZYG')),
            'CHI2': str(record.INFO.get('CHI2')),
            'FISHERPHREDSCORE': str(record.INFO.get('FISHERPHREDSCORE')),
            'INH': str(record.INFO.get('INH')),
            'BESTSTATE': str(record.INFO.get('BESTSTATE')),
            'COVSTATE': str(record.INFO.get('COVSTATE')),
            'SOMATIC': str(record.INFO.get('SOMATIC')),
            'DENOVO': str(record.INFO.get('DENOVO')),
            'GTF_DP': str(record.gt_depths[0]),
            'GTF_AD': str(record.gt_alt_depths[0]),
            'MULTIALLELIC': str(record.INFO.get('OLD_MULTIALLELIC')) or None,
            'AAF': str(float(record.gt_alt_depths[0]) / float(record.gt_depths[0]))}

    return info


def parse_platypus_vcf_record(record):

    multi_allelic = record.INFO.get('OLD_MULTIALLELIC') or False

    if multi_allelic:
        tr = record.INFO.get('TR')[0]
    else:
        tr = record.INFO.get('TR')

    if float(record.INFO.get('TC')) < 1:
        aaf = "0"
    else:
        aaf = str(float(tr) / float(record.INFO.get('TC')))

    info = {'DP': str(tr),
            'FR': str(record.INFO.get('FR')),
            'MMLQ': str(record.INFO.get('MMLQ')),
            'TCR': str(record.INFO.get('TCR')),
            'HP': str(record.INFO.get('HP')),
            'WE': str(record.INFO.get('WE')),
            'WS': str(record.INFO.get('WS')),
            'FS': str(record.INFO.get('FS')),
            'TR': str(tr),
            'NF': str(record.INFO.get('NF')),
            'TCF': str(record.INFO.get('TCF')),
            'NR': str(record.INFO.get('NR')),
            'TC': str(record.INFO.get('TC')),
            'END': str(record.INFO.get('END')),
            'MGOF': str(record.INFO.get('MGOF')),
            'SbPval': str(record.INFO.get('SbPval')),
            'START': str(record.INFO.get('START')),
            'ReadPosRankSum': str(record.INFO.get('ReadPosRankSum')),
            'MQ': str(record.INFO.get('MQ')),
            'QD': str(record.INFO.get('QD')),
            'SC': str(record.INFO.get('SC')),
            'BRF': str(record.INFO.get('BRF')),
            'HapScore': str(record.INFO.get('HapScore')),
            'FILTER': str(record.FILTER),
            'MULTIALLELIC': str(record.INFO.get('OLD_MULTIALLELIC')) or None,
            'AAF': aaf}

    return info


def parse_pindel_vcf_record(record):
    info = {'DP': str(record.gt_depths[0]),
            'END': str(record.INFO.get('END')),
            'HOMLEN': str(record.INFO.get('HOMLEN')),
            'HOMSEQ': str(record.INFO.get('HOMSEQ')),
            'SVLEN': str(record.INFO.get('SVLEN')),
            'SVTYPE': str(record.INFO.get('SVTYPE')),
            'NTLEN': str(record.INFO.get('NTLEN')),
            'GTF_DP': str(record.gt_depths[0]),
            'GTF_AD': str(record.gt_alt_depths[0]),
            'AAF': str(float(record.gt_alt_depths[0]) / float(record.gt_depths[0])),
            'FILTER': str(record.FILTER),
            'MULTIALLELIC': str(record.INFO.get('OLD_MULTIALLELIC')) or None}

    return info


def var_is_rare(variant_data, threshold):
    """Check if variant is rare, as defined by the passed cutoff
    :param variant_data: A GeminiRow for a single variant.
    :type variant_data: GeminiRow.
    :param threshold: Allele frequency rarity threshold.
    :type threshold: float.
    :returns:  bool -- True or False.
    """

    if variant_data.INFO.get('in_esp') != 0 or variant_data.INFO.get('in_1kg') != 0 or variant_data.INFO.get('in_exac') != 0:
        if variant_data.INFO.get('max_aaf_all') > threshold:
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

    if variant_data.INFO.get('cosmic_ids') is not None:
        return True
    else:
        return False


def var_is_in_clinvar(variant_data):
    """Check if variant is in the ClinVar database
    :param variant_data: A GeminiRow for a single variant.
    :type variant_data: GeminiRow.
    :returns:  bool -- True or False.
    """

    if variant_data.INFO.get('clinvar_sig') is not None:
        return True
    else:
        return False


def var_is_pathogenic(variant_data):
    """Check if variant is listed as pathogenic in ClinVar
    :param variant_data: A GeminiRow for a single variant.
    :type variant_data: GeminiRow.
    :returns:  bool -- True or False.
    """

    if variant_data.INFO.get('clinvar_sig') is not None:
        if "pathogenic" in variant_data.INFO.get('clinvar_sig'):
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
    if variant_data.INFO.get('impact_severity') != "LOW":
        return True
    else:
        return False


def var_in_gene(variant_data, genes):
    """Check if variant has a gene name associated with it
    :param variant_data: A GeminiRow for a single variant.
    :type variant_data: GeminiRow.
    :returns:  bool -- True or False.
    """
    if variant_data.INFO.get('gene') in genes:
        return True
    else:
        return False


def var_is_lof(variant_data):
    if variant_data.INFO.get('is_lof'):
        return True
    else:
        return False


def var_is_coding(variant_data):
    if variant_data.INFO.get('is_coding'):
        return True
    else:
        return False


def var_is_splicing(variant_data):
    if variant_data.INFO.get('is_splicing'):
        return True
    else:
        return False


def parse_rs_ids(variant_data):
    if variant_data.INFO.get('rs_ids') is not None:
        return variant_data.INFO.get('rs_ids').split(',')
    else:
        return []


def parse_cosmic_ids(variant_data):
    if variant_data.INFO.get('cosmic_ids') is not None:
        return variant_data.INFO.get('cosmic_ids').split(',')
    else:
        return []
