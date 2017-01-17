"""Backends for shapeshifting RseQC output to pandas dataframes."""
import re
from odo import resource
import pandas as pd


class ParsingError(Exception):
    pass


@resource.register('.*infer_experiment\.txt', priority=20)
def resource_infer_experiment(uri, **kwargs):
    """odo resource to parse infer_experiment results.

    Note
    ----
    The infer experiment file must be named `*infer_experiment.txt`.

    """
    header = [
        'Fraction of reads failed to assign',
        'Fraction of reads on same strand as gene',
        'Fraction of reads on opposite strand as gene',
        'Read type',
        ]

    with open(uri, 'r') as fh:
        data = fh.read().strip()
        try:
            # SingleEnd
            se_regex = (
                'This is SingleEnd Data\n'
                'Fraction of reads failed to determine: (\d\.\d+)\n'
                'Fraction of reads explained by \"\+\+\,\-\-\"\: (\d\.\d+)\n'
                'Fraction of reads explained by \"\+\-\,\-\+\"\: (\d\.\d+)'
                )

            m = re.match(se_regex, data, flags=re.DOTALL).groups()
            values = [float(x) for x in m]
            values.append('SingleEnd')
            return pd.Series(data=values, index=header)
        except:
            pass

        try:
            # PairEnd
            pe_regex = (
                'This is PairEnd Data\n'
                'Fraction of reads failed to determine: (\d\.\d+)\n'
                'Fraction of reads explained by \"1\+\+\,1\-\-,2\+\-\,2\-\+\"\: (\d\.\d+)\n'
                'Fraction of reads explained by \"1\+\-\,1\-\+,2\+\+,2\-\-\"\: (\d\.\d+)'
                )

            m = re.match(pe_regex, data, flags=re.DOTALL).groups()
            values = [float(x) for x in m]
            values.append('PairEnd')
            return pd.Series(data=values, index=header)
        except:
            pass

        raise ParsingError('File did not conform to SE or PE read pattern.')

@resource.register('.*geneBody_coverage\.txt', priority=20)
def resource_geneBody_coverage(uri, **kwargs):
    """odo resource to parse geneBody_coverage results.

    Note
    ----
    The gene body coverage file must be named `*geneBody_coverage.txt`.

    """
    df = pd.read_csv(uri, sep='\t', index_col=0)
    df.index = [x.replace('.prealn', '') for x in df.index]
    df.index.name = 'sample_id'
    return df

@resource.register('.*tin\.txt', priority=20)
def resource_tin(uri, **kwargs):
    """odo resource to parse tin results.

    Note
    ----
    The gene body coverage file must be named `*tin.txt`.
    """
    df = pd.read_csv(uri, sep='\t', index_col=0)
    df.index = [x.replace('.prealn.bam', '') for x in df.index]
    df.index.name = 'sample_id'
    return df

@resource.register('.*bam_stat\.txt', priority=20)
def resource_tin(uri, **kwargs):
    """odo resource to parse bam_stat results.

    Note
    ----
    The gene body coverage file must be named `*bam_stat.txt`.
    """

    header = ['Total_records', 'QC_failed', 'Optical_PCR_duplicate',
              'Non_primary_hits', 'Unmapped_reads', 'non-unique',
              'unique', 'Read_1', 'Read_2', 'Reads_map_plus', 'Reads_map_minus', 'NonSplice_reads',
              'Splice_reads', 'Reads_mapped_in_proper_pairs',
              'ProperPaired_reads_map_to_different_chrom']

    with open(uri, 'r') as fh:
        values = []
        for row in fh:
            row = row.strip()
            try:
                value = re.match(r'.*[\s\:](\d+)', row).groups()[0]
                values.append(int(value))
            except:
                pass

        return pd.Series(data=values, index=header)
