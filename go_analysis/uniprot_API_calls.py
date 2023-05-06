#uniprot API, to download GO terms
import urllib.parse
import urllib.request
import json
import pandas as pd

def uniprot_api(uniprot):
    # Use the uniprot as input to fetch the protein name
    url_template_uniprot = "https://rest.uniprot.org/uniprotkb/{}.json"
    url_uniprot = url_template_uniprot.format(uniprot)

    with urllib.request.urlopen(url_uniprot) as link:
        print(f'Extracting protein name for {uniprot}...')
        data_uniprot = json.loads(link.read().decode())
        #print(data_uniprot)

        try:
            go_list = []
            for count,value in enumerate(data_uniprot['uniProtKBCrossReferences']):
                if data_uniprot['uniProtKBCrossReferences'][count]["database"] == "GO":
                    go_id = data_uniprot['uniProtKBCrossReferences'][count]["id"]
                    #print(go_id)
                    go_list.append(go_id)
        except:
            go_id = "unknown"
            go_list.append(go_id)

        try:
            pfam_list = []
            for count,value in enumerate(data_uniprot['uniProtKBCrossReferences']):
                if data_uniprot['uniProtKBCrossReferences'][count]["database"] == "Pfam":
                    pfam_id = data_uniprot['uniProtKBCrossReferences'][count]["id"]
                    #print(pfam_id)
                    pfam_list.append(pfam_id)
        except:
            pfam_id = "unknown"
            pfam_list.append(pfam_id)


        if len(go_list) == 0:
            go_id = "unknown"
            go_list.append(go_id)
        if len(pfam_list) == 0:
            pfam_id = "unknown"
            pfam_list.append(pfam_id)

        return go_list
        #return go_list, pfam_list

#fetching a list of all IDPs from local machine
idr_all = pd.read_csv('../md_simulations/data/idr_all.csv.gz',header=0,sep=';')
uniprots=idr_all["uniprot"].to_list()
uniprots = sorted(uniprots)
uniprots = sorted(list(set(uniprots)))
#print(len(uniprots))

#making a dictionary with keys = uniprot name, values = assigned GO terms
#this takes a while
uniprot_dict = {}

for uniprot in sorted(uniprots):
    go = uniprot_api(uniprot)
    uniprot_dict[uniprot] = go #was [go, pfam] also above
