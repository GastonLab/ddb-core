import sys
from cyvcf2 import VCF


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
    info = {'FILTER': str(record.FILTER),
            'GTF_DP': str(record.gt_depths[0]),
            'GTF_AD': str(record.gt_alt_depths[0])}

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
            'GTF_AD': str(record.gt_alt_depths[0])}

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
            'GTF_DP': str(record.gt_depths[0])}

    return info


def parse_scalpel_vcf_record(record):
    info = {'AVGCOV': str(record.INFO.get('AVGCOV')),
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
            'GTF_AD': str(record.gt_alt_depths[0])}

    return info