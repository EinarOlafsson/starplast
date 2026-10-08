from pathlib import Path
import sys
sys.path.insert(0,str(Path.cwd()))
from scripts.notebook_runner import ExecutedNotebook
nb=ExecutedNotebook('Host installed identifiers: inspect failed replay association')
nb.md('The replay guard refused an installed-to-native identifier association. Inspect the original installed keys before adapting the comparison; never guess their format or silently omit unmatched genes.')
nb.code("import pandas as pd", "from starplast import organisms as O", "from pathlib import Path", "for organism in (O.HUMAN,O.MOUSE):", "    installed=pd.read_parquet(Path('starplast/data')/O.HOST_TABLES[organism])", "    print(organism,installed.columns.tolist())", "    print(installed[['host_id','host_name']].head().to_string(index=False))")
nb.write('results/host_expression_source_review_2026_10_08/comparison_v2/identifier_inspection.ipynb')
