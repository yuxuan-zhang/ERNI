import numpy as np
import pandas as pd
# import periodictable as pt
# from periodictable import constants
from scipy import interpolate
from bs4 import BeautifulSoup
import urllib.request
import os

data = {}
iso_num = []
df = pd.DataFrame()
df_name = pd.DataFrame()

url = ['http://www.nndc.bnl.gov/exfor/servlet/X4sGetTabData?datfile=E4R21934_tdat.dat&emin=1e-05&emax=2e%2B07']

# url = ['http://www.nndc.bnl.gov/exfor/servlet/X4sGetTabData?datfile=E4R21920_tdat.dat&emin=1e-05&emax=2e%2B08',
#        'http://www.nndc.bnl.gov/exfor/servlet/X4sGetTabData?datfile=E4R21921_tdat.dat&emin=1e-05&emax=2e%2B08',
#        'http://www.nndc.bnl.gov/exfor/servlet/X4sGetTabData?datfile=E4R21922_tdat.dat&emin=1e-05&emax=2e%2B08',
#        'http://www.nndc.bnl.gov/exfor/servlet/X4sGetTabData?datfile=E4R21923_tdat.dat&emin=1e-05&emax=1.5e%2B08'
#        ]

# Co
# url = ['http://www.nndc.bnl.gov/exfor/servlet/X4sGetTabData?datfile=E4R21891_tdat.dat&emin=1e-05&emax=2e%2B07']
# Gd
# url = ['http://www.nndc.bnl.gov/exfor/servlet/X4sGetTabData?datfile=E4R21866_tdat.dat&emin=1e-05&emax=2e%2B07',
#        'http://www.nndc.bnl.gov/exfor/servlet/X4sGetTabData?datfile=E4R21867_tdat.dat&emin=1e-05&emax=2e%2B07',
#        'http://www.nndc.bnl.gov/exfor/servlet/X4sGetTabData?datfile=E4R21868_tdat.dat&emin=1e-05&emax=2e%2B07',
#        'http://www.nndc.bnl.gov/exfor/servlet/X4sGetTabData?datfile=E4R21869_tdat.dat&emin=1e-05&emax=2e%2B07',
#        'http://www.nndc.bnl.gov/exfor/servlet/X4sGetTabData?datfile=E4R21870_tdat.dat&emin=1e-05&emax=2e%2B07',
#        'http://www.nndc.bnl.gov/exfor/servlet/X4sGetTabData?datfile=E4R21871_tdat.dat&emin=1e-05&emax=2e%2B07',
#        'http://www.nndc.bnl.gov/exfor/servlet/X4sGetTabData?datfile=E4R21872_tdat.dat&emin=1e-05&emax=2e%2B07'
#        ]
# U
# url = ['http://www.nndc.bnl.gov/exfor/servlet/X4sGetTabData?datfile=E4R1231_tdat.dat&emin=1e-05&emax=2e%2B07',
#        'http://www.nndc.bnl.gov/exfor/servlet/X4sGetTabData?datfile=E4R1232_tdat.dat&emin=1e-05&emax=3e%2B07',
#        'http://www.nndc.bnl.gov/exfor/servlet/X4sGetTabData?datfile=E4R1233_tdat.dat&emin=1e-05&emax=3e%2B07',
#        'http://www.nndc.bnl.gov/exfor/servlet/X4sGetTabData?datfile=E4R1234_tdat.dat&emin=1e-05&emax=3e%2B07'
#        ]

# Cd
# url = ['http://www.nndc.bnl.gov/exfor/servlet/X4sGetTabData?datfile=E4R21355_tdat.dat&emin=1e-05&emax=2e%2B07',
#        'http://www.nndc.bnl.gov/exfor/servlet/X4sGetTabData?datfile=E4R748_tdat.dat&emin=1e-05&emax=2e%2B07',
#        'http://www.nndc.bnl.gov/exfor/servlet/X4sGetTabData?datfile=E4R749_tdat.dat&emin=1e-05&emax=2e%2B07',
#        'http://www.nndc.bnl.gov/exfor/servlet/X4sGetTabData?datfile=E4R750_tdat.dat&emin=1e-05&emax=2e%2B07',
#        'http://www.nndc.bnl.gov/exfor/servlet/X4sGetTabData?datfile=E4R751_tdat.dat&emin=1e-05&emax=2e%2B07',
#        'http://www.nndc.bnl.gov/exfor/servlet/X4sGetTabData?datfile=E4R752_tdat.dat&emin=1e-05&emax=2e%2B07',
#        'http://www.nndc.bnl.gov/exfor/servlet/X4sGetTabData?datfile=E4R755_tdat.dat&emin=1e-05&emax=2e%2B07',
#        'http://www.nndc.bnl.gov/exfor/servlet/X4sGetTabData?datfile=E4R753_tdat.dat&emin=1e-05&emax=2e%2B07'
#        ]


# url = ['http://www.nndc.bnl.gov/exfor/servlet/E4sGetTabSect?SectID=524330&req=747&PenSectID=810925',
#        'http://www.nndc.bnl.gov/exfor/servlet/E4sGetTabSect?SectID=524426&req=748&PenSectID=811021',
#        'http://www.nndc.bnl.gov/exfor/servlet/E4sGetTabSect?SectID=524479&req=749&PenSectID=811074',
#        'http://www.nndc.bnl.gov/exfor/servlet/E4sGetTabSect?SectID=524540&req=750&PenSectID=811135',
#        'http://www.nndc.bnl.gov/exfor/servlet/E4sGetTabSect?SectID=524637&req=751&PenSectID=811232',
#        'http://www.nndc.bnl.gov/exfor/servlet/E4sGetTabSect?SectID=524686&req=752&PenSectID=811281',
#        'http://www.nndc.bnl.gov/exfor/servlet/E4sGetTabSect?SectID=524755&req=753&PenSectID=811350',
#        'http://www.nndc.bnl.gov/exfor/servlet/E4sGetTabSect?SectID=524968&req=755&PenSectID=811563'
#        ]
main_dir = os.path.dirname(os.path.abspath(__file__))

for n in range(len(url)):
    # Read from html with lxml parser
    html_content = urllib.request.urlopen(url[n]).read()
    soup = BeautifulSoup(html_content, 'lxml')
    p = soup.p
    # get text
    text = soup.get_text()  # <class 'str'>
    # split str by '\r\n'
    list1 = text.splitlines()   # list1 = text.split(sep='\r\n') has same result
    # get where data starts
    data_start = [i for i, s in enumerate(list1) if '1E-05' in s]
    # find isotope names and database used
    ele = [s for i, s in enumerate(list1) if 'NUCLEUS' in s]  # filter(lambda s: 'NUCLEUS' in s, list1)  #
    library = [s for i, s in enumerate(list1) if 'LIBRARY' in s]  # filter(lambda s: 'LIBRARY' in s, list1)  #
    ele_s = ele[0].split()
    lib_s = library[0].split()
    # Isotope name
    iso = ele_s[1].split(sep='-')
    element = iso[0]
    number = iso[1]
    # Database used
    database = lib_s[1]
    # Directory chosen
    if database == 'ENDF/B-VIII.b4':
        _database = 'ENDF_VIII'
    else:
        _database = 'ENDF_VII'
    save_path = main_dir + '/data_web/' + _database + '/'
    # Separate lines
    lst = text.split(sep='Lin-Lin\r\n')
    line1 = lst[0].splitlines()
    lst[0] = line1[-1]
    # Split the column interested
    lst1 = [i.split()[0] for i in lst]  # E,eV  scrap from start
    lst2 = [i.split()[-1] for i in lst]  # Sig,b scrap from end
    del lst1[-1]  # drop the non float part
    del lst2[-1]  # drop the non float part
    # Convert char to float
    lst1 = [float(i) for i in lst1]
    lst2 = [float(i) for i in lst2]
    lst1.insert(0, 'E_eV')
    lst2.insert(0, 'Sig_b')
    lst1.insert(0, element)
    lst2.insert(0, iso[1])
    iso_num.append(int(number))
    df = pd.DataFrame(lst2, index=lst1)
    filename = ''.join(ele_s[1])
    df.to_csv(save_path + ele_s[1]+'.csv', header=False)
    print(lst1[-1])  # check the energy range
    print(iso)
    print(ele)

print(df.head())
# print(df.tail())
print(iso_num)


"""

"""
