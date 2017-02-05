#!/usr/bin/env python3
"""
   getmy23andme.py - Retrieve DNA matches information from 23andMe
   Copyright (C) 2015 Giulio Genovese (giulio.genovese@gmail.com)

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

   Written by Giulio Genovese <giulio.genovese@gmail.com>
"""

import sys, argparse, getpass, time, re, json, html.parser, pandas as pd, csv

try:
    import requests
except ImportError:
    sys.stderr.write('You need to install the requests module first\n')
    sys.stderr.write('(run this in your terminal: "python3 -m pip install requests" or "python3 -m pip install --user requests")\n')
    exit(2)

class Session:
    def __init__(self, username, password, verbose, logfile, timeout):
        self.username = username
        self.password = password
        self.verbose = verbose
        self.logfile = logfile
        self.timeout = timeout
        self.s = requests.Session()
        self.login()

    def login(self):
        url = 'https://www.23andme.com/cas/signin/'
        data = { 'username': self.username, 'password': self.password, '__source_node__': 'start', '__form__': 'login'}
        while True:
            if self.verbose:
                self.logfile.write('[' + time.strftime("%Y-%m-%d %H:%M:%S") + ']: Downloading: ' + url + '\n')
            try:
                r = self.s.post(url, data = data, timeout = self.timeout)
            except requests.exceptions.ReadTimeout:
                if self.verbose:
                    self.logfile.write('[' + time.strftime("%Y-%m-%d %H:%M:%S") + ']: Read timed out\n')
                continue
            except requests.exceptions.ConnectionError:
                if self.verbose:
                    self.logfile.write('[' + time.strftime("%Y-%m-%d %H:%M:%S") + ']: Connection aborted\n')
                time.sleep(self.timeout)
                continue
            if self.verbose:
                self.logfile.write('[' + time.strftime("%Y-%m-%d %H:%M:%S") + ']: Status code: ' + str(r.status_code) + '\n')
            cookies = requests.utils.dict_from_cookiejar(self.s.cookies)
            self.cookies = { 'username': cookies['username'], 'b': cookies['b'], 'uuid': cookies['uuid'], 'session': cookies['session'] }
            return

    def get_url(self, url, xhr = False, data = None):
        headers = { 'X-Requested-With': 'XMLHttpRequest' } if xhr else None
        while True:
            if self.verbose:
                self.logfile.write('[' + time.strftime("%Y-%m-%d %H:%M:%S") + ']: Downloading: ' + url + (' ' + str(data) if data else '') + '\n')
            try:
                if data:
                    r = self.s.post(url, cookies = self.cookies, data = data, headers = headers, timeout = self.timeout)
                else:
                    r = self.s.get(url, cookies = self.cookies, headers = headers, timeout = self.timeout)
            except requests.exceptions.ReadTimeout:
                if self.verbose:
                    self.logfile.write('[' + time.strftime("%Y-%m-%d %H:%M:%S") + ']: Read timed out\n')
                continue
            except requests.exceptions.ConnectionError:
                if self.verbose:
                    self.logfile.write('[' + time.strftime("%Y-%m-%d %H:%M:%S") + ']: Connection aborted\n')
                time.sleep(self.timeout)
                continue
            if self.verbose:
                self.logfile.write('[' + time.strftime("%Y-%m-%d %H:%M:%S") + ']: Status code: ' + str(r.status_code) + '\n')
            try:
                r.raise_for_status()
            except requests.exceptions.HTTPError:
                if self.verbose:
                    self.logfile.write('[' + time.strftime("%Y-%m-%d %H:%M:%S") + ']: HTTPError\n')
                time.sleep(self.timeout)
                continue
            text = html.parser.unescape(r.text)
            if self.verbose and xhr:
                self.logfile.write(text + '\n')
            if r.text == '191919':
                self.login()
            else:
                return text

    # this function retrieves the list of profiles from the https://www.23andme.com/you/ page
    # (maybe there is a more direct way to request this list but I could not figure it out)
    def get_profiles(self):
        text = self.get_url('https://www.23andme.com/you/')
        text = html.parser.unescape(re.sub(' *\n *', '', text))

        regexp = re.compile('dataLayer = \[.*?\];')
        res = regexp.search(text)
        line = text[res.span()[0]:res.span()[1]]
        dataLayer = json.loads(line[12:-1])
        ids = [dataLayer[0]['profile_id']]

        regexp = re.compile('<div class=\"(profile-name|user-name)\">.*?</div>')
        res = regexp.search(text)
        line = text[res.span()[0]:res.span()[1]]
        labels = [re.sub('  *', ' ', re.sub('<.*?>', ' ', line)).strip()]

        regexp = re.compile('<li><a id=\"profile_option_' + '[a-z0-9]' * 16 + '\" class=\"profile_option\" href=\"#\">.*?</a></li>')
        for res in regexp.finditer(text):
          line = text[res.span()[0]:res.span()[1]]
          ids.append(line[26:42])
          labels.append(line[76:-9])

        return (dataLayer, { 'people_ids': ids, 'people_labels': labels })
    
    # this function retrieves the JSON variable embedded in the https://www.23andme.com/you/inheritance/ page
    # (maybe there is a more direct way to request the JSON variable but I could not figure it out)
    def get_inheritance(self):
        text = self.get_url('https://www.23andme.com/you/inheritance/')
        text = html.parser.unescape(re.sub(' *\n *', '', text))
    
        regexp = re.compile('var inheritance = new Inheritance\(\'genome_view\', {.*?}\);')
        res = regexp.search(text)
        line = text[res.span()[0]:res.span()[1]]
        return json.loads(line[49:-2])
    
    def get_relfinder(self, uid):
        text = self.get_url('https://www.23andme.com/you/relfinder/fetch/?profile_id=' + uid, True)
        return json.loads(text)
    
    def get_ibdview(self, p1, p2):
        text = self.get_url('https://www.23andme.com/you/ibdview/', True, { 'p1': p1, 'p2': p2 })
        return json.loads(text)

    def get_ancestry_finder(self, uid):
        text = self.get_url('https://www.23andme.com/you/labs/ancestry_finder/lookup/?profile_id_encrypted=' + uid, True)
        return json.loads(text)

    def get_ancestry_finder_csv(self, uid):
        text = self.get_url('https://www.23andme.com/you/labs/ancestry_finder/export/?profile_id_encrypted=' + uid, True)
        return text # pd.read_csv(io.StringIO(text))

    def get_gender(self, uid):
        if uid in ['v$SP1_FATHER_V4', 'v$SP1_SON1_V2', 'v$SP1_SON2_V2', 'v$SP1_MOTHERS_FATHER_V2', 'v$SP1_FATHERS_FATHER_V2', 'v$NA18558', 'v$NA18944', 'v$NA19160']:
            return 'Male'
        if uid in ['v$SP1_MOTHER_V4', 'v$SP1_DAUGHTER_V2', 'v$SP1_MOTHERS_MOTHER_V2']:
            return 'Female'

        text = self.get_url('https://www.23andme.com/user/?profile=' + uid)
        text = html.parser.unescape(re.sub(' *\n *', '', text))

        regexp = re.compile('<p><strong>Sex:</strong>(Female|Male)</p>')
        res = regexp.search(text)
        if res:
            line = text[res.span()[0]:res.span()[1]]
            return line[24:-4]
        else:
            return 'Unknown'

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'Retrieve DNA matches from 23andMe (26 Jun 2016)', add_help = False, usage = 'getmy23andme.py -u <username> -p <password> [options]')
    parser.add_argument('-u', metavar = '<STR>', type = str, help = '23andMe username [prompt]')
    parser.add_argument('-p', metavar = '<STR>', type = str, help = '23andMe password [prompt]')
    parser.add_argument('-v', action = 'store_true', default = False, help = 'whether to use verbose mode [False]')
    parser.add_argument('-t', metavar = '<INT>', type = int, default = 60, help = 'timeout in seconds [60]')
    parser.add_argument('-o', metavar = '<STR>', type = str, help = 'output prefix [account_id]')
    parser.add_argument('-x', action = 'store_true', default = False, help = 'whether to download inheritance and ibdview tables [False]')
    parser.add_argument('-h', metavar = '<FILE>', type = str, help = 'previously downloaded inheritance table file')
    parser.add_argument('-i', metavar = '<FILE>', type = str, help = 'previously downloaded ibdview table file')
    try:        
        parser.add_argument('-l', metavar = '<FILE>', type = argparse.FileType('w', encoding = 'UTF-8'), default = sys.stderr, help = 'output log file [stderr]')
    except TypeError:
        sys.stderr.write('Python >= 3.4 is required to run this script\n')
        sys.stderr.write('(see https://docs.python.org/3/whatsnew/3.4.html#argparse)\n')
        exit(2)

    # extract arguments from the command line
    try:
        parser.error = parser.exit
        args = parser.parse_args()
    except SystemExit:
        parser.print_help()
        exit(2)

    username = args.u if args.u else input("Enter 23andMe username: ")
    password = args.p if args.p else getpass.getpass("Enter 23andMe password: ")

    # initialize a session with 23andMe server
    session = Session(username, password, args.v, args.l, args.t)

    # download list of profiles handled in the account
    dataLayer, profiles = session.get_profiles()
    out = args.o if args.o else dataLayer[0]['account_id']
    pd.DataFrame(profiles).to_csv(out + '.tsv', sep = '\t', na_rep = 'NA', index = False)

    # download list of DNA matches for each profile
    gender = dict()
    for ehid in profiles['people_ids']:
        relfinder = session.get_relfinder(ehid)
        keys = ['share_status', 'desc', 'eiid', 'rel_upper', 'match_id', 'max_grandparents_same_country', 'patside', 'resend_date', 'segs', 'year', 'ehid', 'birth_country_maternal_gma', 'first', 'res', 'rel_alt', 'url', 'birth_country_maternal_gpa', 'first_initial', 'rel_alg', 'last', 'first_name', 'visible', 'tree_url', 'rel_alg_label', 'disc', 'last_initial', 'pct', 'anc', 'locs', 'matside', 'full', 'pat', 'invitation_status', 'rel_label', 'hide_rel', 'favorite', 'new_share_status', 'birth', 'eid', 'mat', 'rel_range', 'rel_lower', 'img', 'self_reported_ashkenazi', 'rel_user', 'updated', 'can_resend', 'sex', 'birth_country_paternal_gma', 'last_name', 'added', 'surs', 'birth_country_paternal_gpa', 'intro_status']
        df = pd.DataFrame(columns = keys)
        i = 0
        for match in relfinder['matches']:
            if match['ehid']:
                gender[match['ehid']] = match['sex']
            for key, value in match.items():
                df.set_value(i, key, value)
            i += 1
        df.to_csv(out + '.' + ehid + '.relfinder.tsv', sep = '\t', na_rep = 'NA', index = False)

    # download list of IBD segments for each profile (discontinued as of June 2016)
    if False:
      for ehid in profiles['people_ids']:
          ancfinder = session.get_ancestry_finder(ehid)
          keys = ['num', 'ehid', 'label', 'sl', 'nd', 'cm', 'st', 'pgf_ashk', 'pgf_bc', 'pgm_ashk', 'pgm_bc', 'mgf_ashk', 'mgf_bc', 'mgm_ashk', 'mgm_bc']
          df = pd.DataFrame(columns = keys)
          i = 0
          for num, match in ancfinder['ancfinder_result']['af_record'].items():
              for segment in match['segments']:
                  df.set_value(i, 'num', num)
                  if 'public_data' in match:
                      df.set_value(i, 'ehid', match['public_data']['profile_url'][15:])
                      df.set_value(i, 'label', match['public_data']['label'])
                  for key, value in segment.items():
                      df.set_value(i, key, value)
                  for x in 'p', 'm':
                      for y in 'f', 'm':
                          for z in 'ashk', 'bc':
                              df.set_value(i, x + 'g' + y + '_' + z, match['ancestry'][x + 'g' + y][z])
                  i += 1
          df.to_csv(out + '.' + ehid + '.ancfinder.tsv', sep = '\t', na_rep = 'NA', index = False)

    # download match details for each pair of shared profiles
    if args.x:
        if args.h and args.i:
            inheritance = pd.read_csv(args.h, sep = '\t')
            people_ids = set(inheritance['people_ids'])
            df = pd.read_csv(args.i, sep = '\t')
        else:
            keys = ['p2', 'unassayble_regions', 'function_call', 'p1', 'intervals']
            df = pd.DataFrame(columns = keys)
        inheritance = session.get_inheritance()
        inheritance['gender'] = [gender[uid] if uid in gender else session.get_gender(uid) for uid in inheritance['people_ids']]
        pd.DataFrame({ 'people_ids': inheritance['people_ids'], 'people_labels': inheritance['people_labels'], 'gender': inheritance['gender'] }).to_csv(out + '.inheritance.tsv', sep = '\t', na_rep = 'NA', index = False)
        null = '{"20": [[], []], "21": [[], []], "22": [[], []], "1": [[], []], "3": [[], []], "2": [[], []], "5": [[], []], "4": [[], []], "7": [[], []], "6": [[], []], "9": [[], []], "8": [[], []], "Y": [[], []], "X": [[], []], "11": [[], []], "10": [[], []], "13": [[], []], "12": [[], []], "15": [[], []], "14": [[], []], "17": [[], []], "16": [[], []], "19": [[], []], "18": [[], []]}'
        for i in range(len(inheritance['people_ids'])):
            for j in range(0,i):
                p1 = inheritance['people_ids'][j]
                p2 = inheritance['people_ids'][i]
                if not (args.h and args.i and p1 in people_ids and p2 in people_ids):
                    ibdview = session.get_ibdview(p1, p2)
                    if ibdview['intervals'] != null:
                        for key, value in ibdview.items():
                            df.set_value(hash((p1, p2)), key, value)
        df.to_csv(out + '.ibdview.tsv', sep = '\t', na_rep = 'NA', index = False, quoting = csv.QUOTE_NONE)
