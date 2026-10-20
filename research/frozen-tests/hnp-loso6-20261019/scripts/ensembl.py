"""bulk sequence fetch from the Ensembl REST service (GRCh37 or GRCh38), batches of 50 regions."""
import json, time, urllib.request

def fetch(regions, assembly=37):
    host = 'grch37.rest.ensembl.org' if assembly == 37 else 'rest.ensembl.org'
    out = []
    for i in range(0, len(regions), 50):
        req = urllib.request.Request(f'https://{host}/sequence/region/human', data=json.dumps({'regions': regions[i:i + 50]}).encode(),
                                     headers={'Content-Type': 'application/json', 'Accept': 'application/json'})
        for a in range(6):
            try: res = json.load(urllib.request.urlopen(req, timeout=90)); break
            except Exception as e: print('retry', i, a, e, flush=True); time.sleep(2 ** a)
        else: raise SystemExit(f'ensembl batch failed at {i}')
        out += [r['seq'].upper() for r in res]; time.sleep(0.15)
    return out
