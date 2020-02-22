"""AAAAAAAAAAAAAAAAAAAAAAAAAA"""

import re
import sys
import pandas
import requests
from bs4 import BeautifulSoup as soup

def get_opts(subargs: list) -> list:
    "Get all args until end or another parameter"
    opts = []
    for i in subargs:
        if '-' in i:
            break
        else:
            opts.append(i)
    return opts

def get_task(args: list) -> tuple:
    "Get task options from arguments"
    if '-l' in args:
        opts = get_opts(args[args.index('-l')+1:])
        return 'l', opts
    if '-p' in args:
        opts = get_opts(args[args.index('-p')+1:])
        return 'p', opts

def analyze_lines(lines: list, db):
    "Master function for determining mutational similarities accross lines"
    lines, t_lines = process_lines(lines, db)
    if not lines:
        return
    all_muts = get_muts(t_lines, db)
    muts = join_muts(all_muts)
    if not muts:
        expand_search(lines, t_lines, all_muts, db)
    else:
        offer_muts(lines, t_lines, muts, db)

def process_lines(lines: list, db):
    "Filter line names and get full versions from database"
    all_lines = set(db['Tumor'])
    t_lines = []
    d_lines = []
    for l in lines:
        count = 0
        for t in all_lines:
            if l.upper() in t:
                t_lines.append(t)
                break
            else:
                count += 1
        if count == len(all_lines):
            d_lines.append(l)
            lines.pop(l)
    if d_lines:
        print('Lines ', end = "")
        for d in d_lines:
            print(d, end = " ")
        print('not found.')
    return lines, t_lines

def get_muts(t_lines: list, db):
    muts = []
    for t in t_lines:
        muts.append(list(db[db['Tumor'] == t]['Gene']))
    return muts

def join_muts(all_muts):
    set_muts = set(all_muts[0]).intersection(set(all_muts[1]))
    for m in all_muts[2:]:
        set_muts = set_muts.intersection(set(m))
    return list(set_muts)

def offer_muts(lines: list, t_lines: list, muts: list, db):
    "Lets analyze discovered mutations more deeply"
    print('Found common mutations in ', end = '')
    for m in muts:
        print(m, end = ' ')
    print('\n')

    for i in muts:
        mdf = pandas.DataFrame()
        for l in t_lines:
            mdf = mdf.append(db[(db['Gene'] == i) & (db['Tumor'] == l)])
        for d in ['Gene', 'Type', 'OChange', 'GChange']:
            mdf = mdf.drop(d, axis = 1)
        pchanges = mdf['PChange']
        prot = extract_info(request_prot(i))
        mdf['Loc'] = prot.verify_loc(pchanges)
        mdf.set_index('Tumor', inplace = True, drop = True)
        print(prot.name)
        print(prot.descr + '\n')
        print(mdf)

def request_prot(gene):
    r = requests.get(f'https://uniprot.org/uniprot/?query=gene:{gene}+organism:9606&format=txt')
    return r.text

class Prot:
    def __init__(self, name, descr, domains):
        self.name = name
        self.descr = descr
        self.domains = domains
    def verify_loc(self, pchanges):
        loc = []
        for i in pchanges:
            try:
                i = int(filters['ploc'].search(i).group(1))
                diter = iter(self.domains)
                while True:
                    try:
                        d = next(diter)
                        if i in d:
                            loc.append(f'{d.name}: {i - d.start + 1}')
                            del diter
                            break
                    except StopIteration:
                        del diter
                        loc.append('NID')
                        break
            except AttributeError:
                loc.append('-')
            except TypeError:
                loc.append('-')
        return loc

class Domain:
    def __init__(self, start, end, name):
        self.start = start
        self.end = end
        self.name = name
    def __contains__(self, pos):
        if pos >= self.start and pos <= self.end:
            return True
        else:
            return False

filters = {'id': re.compile(r'(?<=RecName: Full=)(.*?)(?=;)'),
           'func': re.compile(r'(?<=FUNCTION: )(.*?)(?= \{)'),
           'span': re.compile(r'(?<=)([0-9]*)(\.\.)([0-9]*)(?=\nFT)'),
           'name': re.compile(r'(?<=/note=")(.*?)(?="\nFT)'),
           'ploc': re.compile(r'[A-Z]([0-9]*)[A-Z]')}

def extract_info(resp):
    resp = resp.replace('\nCC       ', ' ')
    prot = filters['id'].search(resp).group(1)
    descr = filters['func'].search(resp).group(1)
    domains = []
    resp = resp[resp.index('CHAIN'):resp.index('NP_BIND')]
    while True:
        try:
            resp = resp[resp.index('DOMAIN'):]
        except ValueError:
            break
        start, end = filters['span'].search(resp).group(1, 3)
        name = filters['name'].search(resp).group(1)
        domains.append(Domain(int(start), int(end), name))
        resp = resp[10:]
    return Prot(prot, descr, domains)

def expand_search(lines: list, t_lines: list, all_muts: list, db):
    pass

def analyze_prots(prots: list, db):
    "Master function for determining lines with mutations in specified protein(s)"
    prots = process_prots(prots)
    if not prots:
        return

def main():
    mode, opts = get_task(sys.argv[1:])
    db = pandas.read_csv('ccle_optimized.csv')
    if mode == 'l':
        analyze_lines(opts, db)
    elif mode == 'p':
        analyze_prots(opts, db)

if __name__ == '__main__':
    main()
