"""Backends for shapeshifting RseQC output to pandas dataframes."""
import re
from odo import resource
import pandas as pd

def ParsingError(Exception):
    pass

@resource.register('.*infer_experiment\.txt', priority=20)
def resource(uri, **kwargs):
    """odo resource to parse infer_experiment results."""
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
