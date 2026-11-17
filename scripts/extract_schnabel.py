"""Read the source XLSX as OpenXML; never alter its workbook or infer blanks."""
from pathlib import Path
import json,zipfile,xml.etree.ElementTree as ET
ROOT=Path(__file__).resolve().parents[1]
NS={'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
def run():
 p=ROOT/'data/raw/McCreaMorgan/schnabel_Tab3.1.xlsx'
 with zipfile.ZipFile(p) as z:
  strings=[''.join(e.itertext()) for e in ET.fromstring(z.read('xl/sharedStrings.xml')).findall('s:si',NS)]
  out={}
  for i,name in enumerate(('Sheet1','Sheet2','Sheet3'),1):
   rows=[]
   for row in ET.fromstring(z.read(f'xl/worksheets/sheet{i}.xml')).findall('.//s:row',NS):
    values={}
    for cell in row.findall('s:c',NS):
     v=cell.find('s:v',NS)
     if v is None:continue
     col=''.join(c for c in cell.attrib['r'] if c.isalpha());idx=0
     for c in col:idx=idx*26+ord(c)-64
     values[idx-1]=strings[int(v.text)] if cell.attrib.get('t')=='s' else int(float(v.text))
    if values:rows.append([values.get(j) for j in range(20)])
   out[name]=rows
 (ROOT/'data/processed/schnabel_table_raw.json').write_text(json.dumps(out,indent=2)+'\n')
 print('Extracted original cells; blanks preserved.')
if __name__=='__main__':run()
