import json
import requests
import dinopy as dp
from time import sleep

def run_requests(i, res_all, entry_list, c, payload, header, url):
    try:
        for j in range(i, len(entry_list)):
            payload['sequence'] = entry_list[j][0].decode()
            print(payload['sequence'])
            res = requests.post(url, data=json.dumps(payload), headers=header)
            if "errorprob" in res.content.decode():
                res_all.append(res.content.decode())
            c += 1
            if c%1000 == 0:
                print(c)
    except:
        print("Max retries, going to sleep for 100 sec.")
        print("j= " + str(j))
        sleep(100)
        run_requests(j, res_all, entry_list, c, payload, header, url)


if __name__ == '__main__':
    url = 'https://mesa.mosla.de/api/all' #'http://137.248.121.201:5000/api/all'
    header = {'content-type': 'application/json;charset=UTF-8'}
    with open("mesa.json") as json_file:
        config = json.load(json_file)
    f = dp.FastaReader("mcgr_test.fasta")
    payload = config
    payload['asHTML'] = False
    payload["key"] = ''
    c = 0
    res_all = list()
    i = 0
    entry_list = list(f.entries())
    run_requests(i, res_all, entry_list, c, payload, header, url)
    with open("results.txt", "w") as f_:
        for ent in res_all:
            f_.write(ent)
            f_.write("\n")
