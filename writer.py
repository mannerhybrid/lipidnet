from py2neo.database import Graph
from py2neo.data import Node, Relationship
from Bio import Entrez
from bs4 import BeautifulSoup as soup
import datetime
import pickle
import os
Entrez.email = "md.nur.hakim.rosli@gmail.com"
import threading
import time

def worker(retstart):
    retmax=100
    print(threading.currentThread().getName(), "Started!")
    idsseg = get_Idlist(retstart=retstart, retmax=retmax)
    worker_nodes, worker_edges = subgrapher(idsseg, worker_name=threading.currentThread().getName())
    return

def setup():
    author_nodes = []
    pmid_nodes = []
    written_by = []
    
    count = 0
    g = Graph(uri="bolt://localhost:7687", auth=("neo4j", "justjournals"))
    g.schema.create_uniqueness_constraint('PMID','name')
    g.schema.create_uniqueness_constraint('Person','title')
    g.schema.create_uniqueness_constraint('Factor', 'dID')
    g.schema.create_uniqueness_constraint('Qualifier', 'qID')



def get_Idlist(term="Adiposity", retmax=20, retstart="", fullsearch=True):
    Entrez.email = "md.nur.hakim.rosli@gmail.com"

    # 1. eSearch for a list of terms
    idlist = []

    search_handle = Entrez.esearch(db="pubmed",
                                   term=term,
                                   sort="relevance",
                                   retmode="xml",
                                   usehistory="y",
                                   retmax=retmax)
    default_retmax = 100
    default_handle = Entrez.esearch(db="pubmed",term="adiposity", rettype="xml")
    default_soup = soup(default_handle.read(), "xml")
    full_count = int(default_soup.Count.text)
    idlist = [pmid.text for pmid in default_soup.findAll("Id")]
    
    if fullsearch:
        retmax = full_count
        fullsearch = Entrez.esearch(db="pubmed",
                                   term=term,
                                   sort="relevance",
                                   retmode="xml",
                                   usehistory="y",
                                   retmax=retmax)
        fullsoup = soup(fullsearch.read(), "xml")
        idlist = [pmid.text for pmid in fullsoup.find_all("Id")]
        
    print("[+] ID list obtained with {} records.".format(len(idlist)))
    return idlist

def subgrapher(idlist, worker_name="", errant_log=[]):
    import os
    disabled = 0
    passed = 0
    all_nodes = []
    all_edges = []
    id_nodes = []

    mesh_records = []
    author_records = []
    errant_ids = []
    written = 0
    m = open("error_{}.txt".format(worker_name), "w")
    for i in range(len(idlist)):
        fetch_handle = Entrez.efetch("pubmed", id=idlist[i], rettype="xml", retmode="abstract")
        fetch_record = soup(fetch_handle.read(), "xml")
        disabled_score = 0
        try:
            article = fetch_record.PubmedArticleSet.PubmedArticle.MedlineCitation
            meshes = article.findAll("MeshHeading")
            authors = article.Article.AuthorList.findAll("Author")
            year_published = article.DateCompleted.Year.text
            title = article.Article.ArticleTitle.text
            abstract = "\n".join([section.text for section in article.Article.Abstract.findAll("AbstractText")])
            abstractfile = os.path.join(os.getcwd(),
                                        "id_{}.txt".format(idlist[i]))
            thisId = Node("PMID",
                          name="id_{}".format(idlist[i]),
                          id=idlist[i],
                          year=year_published,
                          title="id_{}".format(idlist[i]),
                          full_title=title)
            all_nodes.append(thisId)
            id_nodes.append(thisId)

            try:
                passed += 1
                new_total = len(idlist) - disabled
                author_list = [(idlist[i],
                                str(authors[a].ForeName.text),
                                str(authors[a].LastName.text)) for a in range(len(authors))
                               if authors[a].find("ForeName")]
                qualified_descriptors = [(idlist[i],
                                          meshes[a].QualifierName["UI"],
                                          meshes[a].QualifierName.text,
                                          meshes[a].DescriptorName["UI"],
                                          meshes[a].DescriptorName.text) for a in range(len(meshes))
                                         if meshes[a].find("QualifierName")]
                descriptors = [(idlist[i], "None",
                                "None",
                                meshes[a].DescriptorName["UI"],
                                meshes[a].DescriptorName.text) for a in range(len(meshes))
                               if not meshes[a].find("QualifierName")]
                author_records.extend(author_list)
                mesh_records.extend(descriptors)
                mesh_records.extend(qualified_descriptors)
            except:
                continue

            # print("[+] All ID nodes have been stored!")
            # print(id_nodes)

            full_records = [mesh_records, author_records]
            if worker_name == "":
                with open("data.pickle" 'wb') as f:
                    pickle.dump(full_records, f, pickle.HIGHEST_PROTOCOL)
            else:
                pickle_name = "worker_{}.pickle".format(worker_name)
                with open(pickle_name, 'wb') as f:
                    pickle.dump(full_records, f, pickle.HIGHEST_PROTOCOL)

            # g.schema.create_uniqueness_constraint('Qualifier', 'qID')

            # MeSH Headings
            for mesh_record in mesh_records:
                main_id, qID, qname, dID, dname = mesh_record
                this_Descriptor = Node("Factor",
                                       name=dname.capitalize(),
                                       id=dID,
                                       title=dname.capitalize())
                all_nodes.append(this_Descriptor)
                mainid = "id_{}".format(main_id)
                e = [node['name'] for node in id_nodes].index(mainid)
                node_needed = id_nodes[e]
                descriptor_for = Relationship(this_Descriptor, "IMPORTANT_TO", node_needed)
                all_edges.append(descriptor_for)

                if qID == "None":
                    continue
                else:
                    this_Qualifier = Node("Qualifier",
                                          name=qname.capitalize(),
                                          id=qID,
                                          title=qname.capitalize())
                    all_nodes.append(this_Qualifier)
                    qualified_by = Relationship(this_Qualifier, "QUALIFIER_OF", this_Descriptor)
                    all_edges.append(qualified_by)

            # print("[+] All MeSH nodes and relations have been stored!")

            for author in author_records:
                main_id, firstName, lastName = author
                fullName = firstName + " " + lastName
                fullName = fullName.replace(" ", "_")
                this_Author = Node("Person",
                                   name=fullName,
                                   title=fullName.replace("_", " "))
                mainid = "id_{}".format(main_id)
                e = [node['name'] for node in id_nodes].index(mainid)
                written_By = Relationship(id_nodes[e], "WRITTEN_BY", this_Author)
                all_nodes.append(this_Author)
                all_edges.append(written_By)

            # print("[+] All author nodes and relations have been stored!")

            # Write abstract to file if path doesn't exist
            if not os.path.exists(abstractfile):
                with open(abstractfile, "wb") as w:
                    w.write(title.encode("utf-8"))
                    w.write("\n\n".encode("utf-8"))
                    w.write(abstract.encode("utf-8"))
                    written += 1
                    if i % 5 == 0:
                        print("[+] Worker %s has successfully passed %d out of %d IDs." % (
                        worker_name, written, new_total))

        except:
            disabled_score += 1
            if disabled_score == 0:
                m.write("Errant ID \t Total failed: \t Worker Name")
                m.write(str(idlist[i]) + '\t'+ str(disabled_score) + '\t' + worker_name)
            else:
                m.write(str(idlist[i]) + '\t' + str(disabled_score) + '\t' + worker_name)
            errant_log.append(idlist[i])
            print("[-] {} errant IDs detected in {}.".format(disabled_score, worker_name))
            continue


    m.close()
    return all_nodes, all_edges

def create_subgraph(subgraphs):
    for subgraph in subgraphs:
       try:
           g.create(subgraph)
           print("[+] {} created!".format(subgraph))
       except:
           print("[-] {} already exists!".format(subgraph))

def get_details(ID):
    print(ID)
    fetch_handle = Entrez.efetch("pubmed", id=ID, rettype="xml", retmode="abstract")
    fetch_record = soup(fetch_handle.read(), "xml")
    article = fetch_record.PubmedArticleSet.PubmedArticle.MedlineCitation
    year_published = article.DateRevised.Year.text
    title = article.Article.ArticleTitle.text
    return [title, year_published]

def get_linked_ids(idlist, iterations=1):
    target_id_nodes = []
    linked_id_nodes = []
    all_link_nodes = []
    all_link_edges = []
    idlist2 = []
    for i in range(len(idlist)):
        linker_handle = Entrez.elink(db="pubmed",
                                     dbfrom="pubmed",
                                     id=idlist[i],
                                     rettype="xml")
        linker_record = soup(linker_handle.read(), "xml")
        thisId = Node("PMID",
                      name="id_{}".format(idlist[i]),
                      id=idlist[i],
                      year=get_details(idlist[i])[1],
                      title="id_{}".format(idlist[i]),
                      full_title=get_details(idlist[i])[0]
                      )
        all_link_nodes.append(thisId)

        linker_records = linker_record.findAll("Link")
        linked_id_nodes_now = [Node("PMID",
                                    name="id_{}".format(linker_records[l].Id.text),
                                    id=linker_records[l].Id.text,
                                    year=get_details(linker_records[l].Id.text)[1],
                                    title="id_{}".format(linker_records[l].Id.text),
                                    full_title=get_details(linker_records[l].Id.text)[0]) for l in range(len(linker_records))]

        idlist2.extend(list(set([l.Id.text for l in linker_record.findAll("Link")])))
        all_link_nodes.extend(linked_id_nodes)
        linked_to = [Relationship(thisId, "LINKED_TO", thisNode) for thisNode in linked_id_nodes_now]
        all_link_edges.extend(linked_to)

        if iterations == 1:
            return all_link_nodes, all_link_edges, idlist2
        elif iterations > 1:
            all_link_nodes = []; all_link_edges = []; idlist3 = []
            for j in range(iterations):
                print("[+] Iteration number: {}".format(j))
                nodes, edges, idlist_j = get_linked_ids(idlist)
                all_link_edges.extend(edges); all_link_nodes.extend(nodes)
                idlist3.extend(list(set(idlist_j)))
            return all_link_nodes, all_link_edges, idlist3

def main():
    # setup()
    pwd = os.getcwd()
    current_date = str(datetime.date.today()).replace("-","")
    final_dir = os.path.join(os.getcwd(), current_date)
    if os.path.exists(final_dir):
        os.chdir(final_dir)
    else:
        os.mkdir(final_dir)
        os.chdir(final_dir)
    default_retmax = 1000
    full_count_handle = Entrez.esearch(db="pubmed",term="adiposity", rettype="xml")
    full_count = soup(full_count_handle.read(), "xml")
    full_count = int(full_count.Count.text)
    num_workers = full_count // default_retmax + 1

    workers = []
    retstart = 0
    start = time.time()
    for i in range(num_workers):
        rets = i*default_retmax
        w = threading.Thread(name="Worker_{}".format(i),
                             target=worker,
                             args=(rets,))
        workers.append(w)
        w.start()
    end = time.time()
    print("Completed in {}".format(end - start))

    # idlist = get_Idlist(retmax=1000)
    # nodes, relations = subgrapher(idlist=idlist)
    # create_subgraph(nodes)
    # create_subgraph(relations)


    os.chdir(pwd)

if __name__ == '__main__':
    main()