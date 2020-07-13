from Bio import Entrez
from Bio import Medline
from tqdm import tqdm
import pandas as pd
pd.set_option('display.max_colwidth', -1)
import numpy as np
#Modules for bigram analysis
import nltk
from nltk import bigrams, trigrams
from nltk.corpus import stopwords
from nltk.tokenize import sent_tokenize, word_tokenize 
import collections
import re
#Graphing modules
import matplotlib.pyplot as plt
import networkx as nx
import seaborn as sns
from collections import Counter


# Change this email to your email address
Entrez.email = "jwang@fortressbiotech.com"

people = input('please enter a search: ')

def lenth_results(x):
    Entrez.email=Entrez.email
    keyword = x
    handle = Entrez.esearch(db ='pubmed',
                            retmax=10,
                            retmode ='text',
                            term = keyword)
    results= Entrez.read(handle)
    len_results = results['Count']
    print('Total number of publications that contain the term {}: {}'.format(keyword, results['Count']))    
    return int(len_results)


def search(x):
    Entrez.email=Entrez.email
    pubs = []
    keyword = x
    handle = Entrez.esearch(db ='pubmed',
                            retmax=int(len_results),
                            retmode ='text',
                            mindate = 2005,
                            maxdate = 2020,
                            usehistory='y',
                            idtype = 'acc',
                            term = keyword)
    results= Entrez.read(handle)
    pubs.append(results)
    print('Total number of publications that contain the term {}: {}'.format(keyword, results['Count']))    
    return pubs

def history(x):
    webenv = []
    for x in pubs:
        for i in [x['WebEnv']]:
            webenv.append(i)
    return webenv
def key(x):     
    query_key = []
    for x in pubs:
        for i in [x['QueryKey']]:
            query_key.append(i)
    return query_key

def id_list(x):
    id_list = []
    for x in pubs:
        for i in [x['IdList']]:
            id_list.append(i)
    id_list = id_list[0] #Right now it's a numpy array, turn it into a 1d list
    return id_list

def count(x):     
    count = []
    for x in pubs:
        for i in [x['Count']]:
            count.append(int(i))
    return count

def search_result(history, key, count):
    search_results = {'History':[], 'QueryKey':[], 'Count':[]}
    search_results['History'].append(webenv)
    search_results['QueryKey'].append(query_key)
    search_results['Count'].append(count)
    #Turn to list to more easily iterate over 
    search_results = [search_results]
    search_results2 = pd.DataFrame(pubs)
    lst_col = 'IdList' #Set the column you want to expand 
    search_results2 = pd.DataFrame({
      col:np.repeat(search_results2[col].values, search_results2[lst_col].str.len())
      for col in search_results2.columns.drop(lst_col)}
    ).assign(**{lst_col:np.concatenate(search_results2[lst_col].values)})[search_results2.columns]
    return search_results2


def batch(x):
    ids = list(search_results2['IdList'])
    batch_size = 100
    batches = [ids[x: x+100] for x in range(0, len(ids), batch_size)]
    return batches

def fetch(x):
    record_list = []
    batch_size = 100
    for x in tqdm(batches):
        handle = Entrez.efetch(db = 'pubmed', 
                               rettype = 'medline', 
                               retmode = 'text',
                               id = x,
                               retmax = batch_size,
                               webenv = webenv,
                               query_key = query_key,
                               idtype = 'acc')
        records = Medline.parse(handle)
        record_list.extend(list(records))
        print('Complete')   
    return record_list

    
def citations(x):
    ids = list(search_results2['IdList'])
    batch_size = 100
    batches = [ids[x: x
                   + 100] for x in range(0, len(ids), batch_size)]
    citation_list = []
    for i in tqdm(batches):
        handle = Entrez.elink(dbfrom="pubmed", 
                              db="pmc",
                              LinkName="pubmed_pmc_refs",
                              id=i,
                              retmax = batch_size)
        citations = Entrez.read(handle)
        citation_list.extend(list(citations))
        print('Complete')
    return citation_list


def publication_volume(x):
    #Count number of authors, and weigh by impact factor 
    df = pd.DataFrame(record_list, columns = ['FAU', 'TA', 'AD', 'TI'])
    #This splits up the list into searate rows
    df['FAU'] = df['FAU'].apply(pd.Series) \
        .merge(df, left_index = True, right_index = True
           )
    df['TA'] = df['TA'].str.lower()

    #Create weighting for influential journals
    mask = {'nature':40,
            'lancet':40,
            'jama': 80, 
            'n engl j med': 80,
            'nat commun': 5,
            r'(?s).*':1
            }

    #Basically replicates the journal names column as a new column 
    df['weight'] = df['TA'].replace(mask, regex = True)
    df['weight'] = df['weight'].fillna(1)
    df = df.reindex(df.index.repeat(df.weight))
    authors_flat = [i for i in list(df['FAU'].values.tolist())]

    n = 100
    top_authors = pd.DataFrame.from_records(
        Counter(authors_flat).most_common(n), columns=["Name", "Count"]
        )


    sns.barplot(x = 'Count', y = 'Name', data = top_authors[:20], palette = 'RdBu_r')
    plt.show()
    return plt


        
if __name__ == '__main__':
    len_results = lenth_results(people)
    pubs = search(people) #search
    webenv = history(pubs)#Create a list of WebEnv historys & Query Keys
    query_key = key(pubs)
    id_list = id_list(pubs)
    count = count(pubs) 
    search_results2 = search_result(webenv, query_key, count) #separate out query key/webenv to iterate over for parsing and downloading results
    batches = batch(search_results2) #break into batches of 100
    record_list = fetch(batches) #fetch the records
    data_set = pd.DataFrame(record_list)   #put into a dataframe
    publication_volume(record_list)   #graph out the authors by publication volume
    
    citation_list = citations(batches) #get the citations using the same batches
    citations = pd.DataFrame(citation_list, columns = ['LinkSetDb'])
    citation2 = citations['LinkSetDb'].to_list()
    citation3 = pd.DataFrame(citation2, columns = ['Link'])
    res = pd.concat([citation3,citation3.Link.apply(pd.Series)],axis=1) #separating out the ids nexted within the dictionary link
    res.columns = ['original', 'list', 'dbname', 'pub'] #need to rename columns since some of them are the same now
    res['list']= res['list'].replace(np.nan, str('a')) #change nans to strings, to differentiate from the lists of citation ids
    citations4 = res['list'].to_list()
    lst = []
    for x in citations4:
        if type(x) == list:
            lst.append(len(x))
        else:
            lst.append(0)
    
    senior_authors = data_set['FAU'].to_list()
    sr_au = []
    for x in senior_authors:
        if type(x) == list:
            sr_au.append(x[-1:])
        else:
            sr_au.append(str(x))
    sr_df = pd.DataFrame(sr_au, columns = ['Authors', 'x', 'y'])#Turn list of lists
    sr_au = [i for i in sr_df['Authors']] #..into a list
    
    final_data = pd.DataFrame(sr_au, columns = ['Authors'])
    final_data['publications'] = data_set['PMID']
    final_data['citations'] = lst
    publication_cnts = []
    for x in data_set['PMID']:
        publication_cnts.append(1)
    final_data['publication counts'] = publication_cnts
    
    '''
    TESTING
    '''
    asdf = final_data.groupby('Authors')['citations'].sum(inplace = True).reset_index()
    asdf2 = final_data.groupby('Authors')['publication counts'].sum(inplace = True).reset_index()
    asdf['publication volume'] = asdf2['publication counts']
    #Since the data is still in the same order as when you parsed it, append
    #the data about institutions/contact info so you can extract the emails
    asdf['emails'] = data_set['AD']
    asdf['emails'] = asdf['emails'].apply(lambda x: re.findall(r'[\w\.-]+@[\w\.-]+', str(x)))
    asdf = asdf.sort_values(by = 'citations', ascending = False)
    plt.figure(figsize = (16,8))
    sns.barplot(x = 'citations', y = 'Authors', data = asdf[:20], palette='coolwarm')
    print(asdf[['Authors', 'emails']][:10])
    asdf.to_excel('C:\Users\jwang\Desktop\Python\Biopython\CSV\test2.xlsx')
    
    
# use the scatter function
plt.figure(figsize=(20,14))
ax = sns.scatterplot(x = 'publication volume', y = 'citations', data = asdf, 
                     size = 'publication volume', sizes = (40,400), alpha = 0.5, hue = 'citations', palette = 'muted', legend = False)
def label_point(x, y, val, ax):
    a = pd.concat({'x': x, 'y': y, 'val': val}, axis=1)
    for i, point in a.iterrows():
        ax.text(point['x']+.02, point['y'], str(point['val']))
label_point(asdf['publication volume'], asdf['citations'], asdf['Authors'], ax)


#Show publications over time
top_authors = list(asdf['Authors'].loc[:20])
time_set = pd.DataFrame(record_list, columns = ['FAU', 'EDAT'])
time_set['Year'] = time_set['EDAT'].astype(str).str[0:4].astype(int)
senior_authors = time_set['FAU'].to_list()
sr_au = []
for x in senior_authors:
    if type(x) == list:
        sr_au.append(x[-1:])
    else:
        sr_au.append(str(x))
time_set['Names'] = sr_au
time_set['Names'] = [','.join(map(str, l)) for l in time_set['Names']] #convert list of author names to string
time_set[time_set['Names'].isin(top_authors)] #filter for only authors with lots of citations
time_set = time_set.drop(columns = ['EDAT', 'FAU'])
time_set['count'] = 1
time_set[:10].pivot_table('count', 'Year', 'Names', aggfunc = 'count').plot(
    kind = 'line', marker = 'o', xticks = time_set.Year.unique(), legend = False)


count = time_set.groupby(['Year', 'Names'])['count'].size().reset_index() #Group by years, and count the occurance of number of publications
plt.figure(figsize = (30,20))
df1 = pd.pivot_table(count, values='count', index='Year', columns='Names')
ax = df1.plot(kind='line', legend = False)
label_point(count['count'], count['Names'], count['Year'], ax)

#Number of publications on the subject over the years
data_set.dropna(subset=['EDAT'], inplace=True)
data_set["Year"] = (
    data_set["EDAT"].astype(str).str[0:4].astype(int)
)
yearly = pd.DataFrame(data_set["Year"].value_counts().reset_index())
yearly.columns = ["Year", "Count"]
pal = sns.color_palette("Greens_d", len(yearly))
plt.figure(figsize = (16, 8))
sns.barplot(x="Year", y="Count", data=yearly, palette = np.array(pal[::-1]))
plt.title("Publications over Time")
plt.show()


'''
filter for only key terms you care about
'''


title_table= pd.DataFrame(record_list, columns = ['TI', 'AB', 'AD']) #Pull out title, abstract, contact info
clean_table = title_table
clean_table['AB'] = clean_table['AB'].str.lower()
clean_table['TI'] = clean_table['TI'].str.lower()
search_terms = ['phase', 'trial', 'patients']
#Since it's a dataframe, we need to convert a list of strings into a series that can be applied to the dataframe
clean_table= clean_table[clean_table['TI'].str.contains('|'.join(search_terms))] #Make sure it's a clincial study
#clean_table = clean_table[~clean_table.TI.str.contains('Newly, newly')]  #Turn this on if you want to filter out studies looking at front line therapy patients
#obv_study_terms = ['demographic','association', 'associated', 'retrospective', 'prognostic', 'prediction', 'detection']  
#clean_table = clean_table[~clean_table['AB'].str.contains('|'.join(obv_study_terms))]
clean_table['AB'] = clean_table['AB'].apply(lambda x: re.findall(r'results:.*|findings:.*', str(x))) # remove all the text ahead of results
clean_table = clean_table[clean_table.astype(str)['AB'] != '[]'] #Remove all the blank cels that result 
#clean_table = clean_table[clean_table['AD'].str.contains('@')] #Keep only results that have an email address
clean_table['AD'] = clean_table['AD'].apply(lambda x: re.findall(r'[\w\.-]+@[\w\.-]+', str(x))) #Remove all institution/affiliation info
#Save as a CSV so that if you need to pull up the data again, you don't have to pull from pubmed during working hrs
#Extract results number after results text
#turn abbreviation into full word for next line
df = clean_table

#convert column to a string
df['AB'] = df['AB'].astype(str)
df['TI'] = df['TI'].astype(str)
df['drug names'] = df['TI'].str.strip('[]') #remove the brackets, so it's recognized as a string
df['drug names'] = df['drug names'].str.replace(r'NCT\d+$', '')
df['drug names'] = df['drug names'].str.extract(r'([A-Za-z]{2,6}\d{2,8}\b|[A-Za-z]{2,6}-\d{2,8}\b|\w+mab\b | \w+nib\b |\w+nib\b|\w+sertib\b|\w+stat\b|\w+staurin\b|\w+bine\b|\w+cef\b|\w+asone\w+bital\b|\w+cycline\b|\w+zole\b|\w+platin\b|\w+cept\b|\w+tedin\b)')
#df['drug names'] = df['drug names'].str.extract(r'([^nct][a-z]{1,4}\d{2,7})|\w+stat\b|')
df['trial ID'] = df['AB']
df['trial ID'] = df['trial ID'].str.extract(r'(nct\d+)')

df.to_excel('C:/Users/jwang/Desktop/Python/Biopython/CSV/test3.xlsx')
'''
Parse out late stage drug resuls
'''

df['mOS'] = df['AB'].str.replace('os', 'overall survival')
#df['mOS'] = df['AB'].str.extract('(?<=overall)*[ a-zA-Z]* (\d+.\d+ month)')
#df['mOS'] = df['AB'].str.extract('(?<=overall survival).*? (\d+.\d+ month)')
#df['mOS'] = df['AB'].str.extract('(?<=overall survival).*? (\d+.\d+ month)')
df['mOS'] = df['mOS'].str.extract('(?<=overall)*[^/]* (\d+.\d+ month)')
df['mOS']= df['mOS'].str.strip('month')


df['PFS'] = df['AB'].str.replace('pfs', 'progression free survival')
#df['PFS'] = df['AB'].str.extract('(?<=progression free survival)*[^/]* (\d+.\d+ month)')
df['PFS'] = df['PFS'].str.extract('(?<=progression free).*? (\d+.\d+ month)')
df['PFS']= df['PFS'].str.strip('month')

'''
Parse out early stage drug results
'''
df['ORR'] = df['AB'].str.replace('orr | overall response', 'objective response')
#df['ORR'] = df['AB'].str.extract('(?<=objective response).+?(\d.+?%)')
#df['ORR'] = df['AB'].str.extract('(?<=objective response).+?(\d.|$)%')
df['ORR'] = df['ORR'].str.extract('(?<=objective response).*?(\d.|$%)')
df['ORR']= df['ORR'].str.strip('%')
ORR_Results = df[['TI', 'AB', 'AD', 'drug names', 'ORR']].sort_values('ORR', ascending = False)
ORR_Results.to_excel('C:/Users/jwang/Desktop/Python/Biopython/CSV/ORR_Results.xlsx')

df['PR'] = df['AB'].str.replace('partial response', 'pr')
df['PR']= df[~df['PR'].str.contains('pruritus', na=False)]
#df['ORR'] = df['AB'].str.extract('(?<=objective response).+?(\d.+?%)')
#df['ORR'] = df['AB'].str.extract('(?<=objective response).+?(\d.|$)%')
df['PR'] = df['PR'].str.extract('(?<=pr).*?(\d.|$%)')
df['PR']= df['PR'].str.strip('%')

df['CR'] = df['AB'].str.replace('complete response', 'cr')
#df['ORR'] = df['AB'].str.extract('(?<=objective response).+?(\d.+?%)')
#df['ORR'] = df['AB'].str.extract('(?<=objective response).+?(\d.|$)%')
df['CR'] = df['CR'].str.extract('(?<=cr).*?(\d.|$%)')
df['CR']= df['CR'].str.strip('%')

#Results table for survival
survival = pd.DataFrame()
survival = survival.append(df['drug names'])
survival =survival.append(df['mOS'])
survival = survival.T.dropna()
survival = survival[~survival.mOS.str.contains('-')]
survival = survival.T
print(survival)

progression = pd.DataFrame()
progression = progression.append(df['drug names'])
progression =progression.append(df['PFS'])
progression = progression.T.dropna()
progression = progression[~progression.PFS.str.contains('-')]
progression = progression.T
print(progression)

#Results table for tumor shrinkage
shrink = pd.DataFrame()
shrink = shrink.append(df['drug names'])
shrink = shrink.append(df['ORR'])
shrink = shrink.T.dropna()
shrink = shrink[shrink['ORR'].apply(lambda x: str(x).isdigit())]
shrink = shrink.T
print(shrink)

#Results for partial and complete response
Response = pd.DataFrame()
Response = Response.append(df['drug names'])
Response = Response.append(df['PR'])
Response = Response.append(df['CR'])
Response = Response.fillna(0)
Response = Response.T
Response = Response[Response['PR'].apply(lambda x: str(x).isdigit())]
Response = Response[Response['CR'].apply(lambda x: str(x).isdigit())]
print(Response)


survival.to_excel('C:/Users/jwang/Desktop/Python/Biopython/CSV/survival.xlsx')
progression.to_excel('C:/Users/jwang/Desktop/Python/Biopython/CSV/progression.xlsx')
shrink.to_excel('C:/Users/jwang/Desktop/Python/Biopython/CSV/shrink.xlsx')  
Response.to_excel('C:/Users/jwang/Desktop/Python/Biopython/CSV/response.xlsx')

late_stage = pd.read_excel('C:/Users/jwang/Desktop/Python/Biopython/CSV/survival.xlsx', sheet_name = 'Sheet1', header = 1, index_col =0)
progression = pd.read_excel('C:/Users/jwang/Desktop/Python/Biopython/CSV/progression.xlsx', sheet_name = 'Sheet1', header = 1, index_col = 0)
early_stage = pd.read_excel('C:/Users/jwang/Desktop/Python/Biopython/CSV/shrink.xlsx', sheet_name = 'Sheet1', header = 1, index_col =0)
response = pd.read_excel('C:/Users/jwang/Desktop/Python/Biopython/CSV/response.xlsx', sheet_name = 'Sheet1', header = 0, index_col = 0)

'''
Overall Survival
'''
OS = late_stage.loc['mOS']
#OS = OS.dropna()
OS = OS.sort_values(ascending = True)
OS.plot(kind = 'barh', figsize = (8,6), fontsize = 12, cmap = plt.get_cmap('viridis'), title = 'mOS', legend = False)
#OS.iloc[np.r_[1:10, 15:20].plot(kind = 'barh', figsize = (8,6), fontsize = 12, cmap = plt.get_cmap('viridis'), title = 'mOS', legend = False)
plt.show()
'''
PFS
'''
PFS = progression.loc['PFS']
#PFS = PFS.dropna()
PFS = PFS.sort_values(ascending = True)
#PFS.iloc[np.r_[49:59,85:89]].plot(kind = 'barh', figsize = (15, 10), fontsize = 12, cmap = plt.get_cmap('plasma'), title = 'PFS', legend = False)
PFS.plot(kind = 'barh', figsize = (15, 10), fontsize = 12, cmap = plt.get_cmap('plasma'), title = 'PFS', legend = False)
plt.show()
'''
Objective Response Rate
'''
ORR = early_stage.loc['ORR']
ORR = ORR.dropna()
ORR = ORR.sort_values(ascending = False)
#ORR.iloc[np.r_[1:5, 10:15]].plot(kind = 'bar', figsize = (8, 6),fontsize = 12, cmap = plt.get_cmap('viridis'), title = 'Overall Response Rate', legend = False)
ORR[:15].plot(kind = 'bar', figsize = (8, 6),fontsize = 12, cmap = plt.get_cmap('viridis'), title = 'Overall Response Rate', legend = False)
plt.show()

#######################
'''
CR + PR
'''
CR = response.loc['CR'].fillna(0)
PR = response.loc['PR'].fillna(0)
ORR = pd.concat([CR, PR], axis = 1)
ORR['sum'] = ORR['CR']+ORR['PR']
ORR['sum'] = ORR[ORR['sum'] !='-']
ORR = ORR.sort_values(by = 'sum', ascending = False)
plt.rc('xtick', labelsize=20) 
ORR['sum'].plot(kind = 'bar', stacked = True, figsize = (10,8))
plt.show()
