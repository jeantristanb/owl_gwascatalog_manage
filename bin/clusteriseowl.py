import re
import sys
#p = re.compile("name (.*) is valid")
#>>> result = p.search(s)
#>>> result
#<_sre.SRE_Match object at 0x10555e738>
#>>> result.group(1)
class efo:
  def __init__(self, header,readfile, patternClass, patternSubClass, patternOnproperty, pattern_gwastrait):
         self.list_owlclass=[]
         self.info=[]
         self.dic_ancestor={}
         self.dic_descent={}
         self.sub_class=[]
         self.readf=readfile
         self.patternSubClass=patternSubClass
         self.patternClass=patternClass
         self.patternOnproperty = patternOnproperty
         self.onProperty = []
         self.annotatedTarget = []
         self.pattern_gwastrait= pattern_gwastrait 
         self.annotatedTarget = []
         self.term_replaced_by=None
         res=patternClass.search(header)
         self.is_gwas_trait= "false"
         try :
          self.key=res.group(1)
         except :
            self.key=None 
         self.read_owl()
         self.defined_info()

  def read_owl(self) :
        beginread=False
        for line in self.readf:
            if '<owl:Class' in line :
                 self.list_owlclass.append(efo(line,self.readf, self.patternClass, self.patternSubClass,self.patternOnproperty,self.pattern_gwastrait))
            elif '</owl:Class' in line :
                return self
            else :
              self.info.append(line) 
        print("error not </owl:Class for "+self.key)
        sys.exit(2)
  def extract_re(self, pat,info, baliseneed) :
         subclass=pat.search (info)
         try :
             tmpre=pat.search(info)
             return tmpre.group(1)
         except :
             if baliseneed==True :
               #print(pat, info)
               sys.exit('problem exit')
               exit(2)
             return None
  def extract_pattern(self, cmtinfo, pattern,pattern_re,endbalise):
     balise=True
     list_pattern=[]
     while endbalise not in self.info[cmtinfo]: 
         info=self.info[cmtinfo]
         if pattern in info :
           list_pattern.append(self.extract_re(pattern_re,info, True))
         cmtinfo+=1
     return (cmtinfo, list_pattern)
  def defined_info(self) :
       nbinfo=len(self.info)
       cmtinfo=0
       while cmtinfo < nbinfo :
         info=self.info[cmtinfo]
         if "<rdfs:subClassOf" in info :
           subclass=self.extract_re(self.patternSubClass,info , False)
           if subclass :
             self.sub_class.append(subclass)
           else :
             (cmtinfo,subclass)=self.extract_pattern(cmtinfo, "<owl:onProperty ",self.patternOnproperty,"</rdfs:subClassOf>")
             self.sub_class+=subclass
         elif "<owl:onProperty" in info :
           self.onProperty.append(self.extract_re(self.patternOnproperty,info, True))
         elif "<efo:gwas_trait" in info :
           self.is_gwas_trait=self.extract_re(self.pattern_gwastrait,info, True)
           self.is_gwas_trait.lower()
         elif "<efo:reason_for_obsolescence" in info :
           spl=info.split('>')   
           if len(spl)>1 :
              if 'http' in spl[1] :
                try  :
                  newhttp=re.search("(?P<url>http?://[^\s]+[0-9])", spl[1]).group("url")
                  self.term_replaced_by=newhttp
                except :
                  print('error \n'+info+'\n'+spl[1])
         cmtinfo+=1
  def get_key(self) :
         return self.key
  def get_annotated_target(self) :
      return self.annotatedTarget
  def is_gwastrait (self) :
      if self.is_gwas_trait == 'true' :
          return True
      elif self.is_gwas_trait == 'false' :
          return False
      print(' self.is_gwas_trait '+self.is_gwas_trait+' not know')
      sys.exit(2)
  def add_ancestor(self, ancestor):
       self.dic_ancestor[ancestor.get_key()]=ancestor
  def add_descent(self, descent):
       self.dic_descent[descent.get_key()]=descent
  def return_desc_rec(self) :
      list_descent=[]
      for descent in self.dic_descent.keys() :
           list_descent.append(descent)
           list_descent+=self.dic_descent[descent].return_desc_rec()
              
      return list_descent

class owl :
    def __init__(self, File) :
        self.patternclass = re.compile("<owl:Class rdf:about=\"(.*)\"")
        #        <rdfs:subClassOf rdf:resource="http://www.ebi.ac.uk/efo/EFO_0000508"/>
        self.patternSubClass= re.compile("<rdfs:subClassOf rdf:resource=\"(.*)\"")
        self.obscolescence={}
        self.patternOnproperty = re.compile("<owl:onProperty rdf:resource=\"(.*)\"/>")
        #        <efo:gwas_trait rdf:datatype="http://www.w3.org/2001/XMLSchema#boolean">true</efo:gwas_trait>
        #        <efo:gwas_trait rdf:datatype="http://www.w3.org/2001/XMLSchema#boolean">true</efo:gwas_trait>
        self.pattern_gwastrait= re.compile("<efo:gwas_trait.*>(.*)</efo:gwas_trait>")
        self.readf=open(File)
        self.list_efo=[]
        self.dic_all={}
        self.end_file=False
        self.read_all()
    def read_next_efo(self):
        if not self.readf :
            self.end_file=True
            return None
        new_owl=[]
        for line in self.readf :
            if '<owl:Class' in line :
                efo_res=efo(line, self.readf, self.patternclass,self.patternSubClass, self.patternOnproperty,self.pattern_gwastrait)
                self.list_efo.append(efo_res)
                if efo_res.get_key() :
                   self.dic_all[efo_res.get_key()]=efo_res
                new_term=efo_res.term_replaced_by
                if new_term :
                    self.obscolescence[new_term]=efo_res.get_key()
                return efo_res
    def read_all(self):      
       balisegood=True
       while balisegood :
         efo_val=self.read_next_efo()
         if not efo_val:
             return()
         #print(efo_val.get_key(),efo_val.is_gwastrait())
    def get_allkey(self) :
        list_efo=[]
        for efol in self.list_efo :
            key=efol.get_key()
            if key :
                list_efo.append(key)
        return list_efo
    def build_tree(self) :
        for key in self.dic_all.keys() :
            efo_val=self.dic_all[key]
            if len(efo_val.sub_class)>0 :
               for x in efo_val.sub_class :
                   if x in self.dic_all:
                       efo_val.add_ancestor(self.dic_all[x]) 
                       self.dic_all[x].add_descent(efo_val)
                   else :
                        print 'not found '+x+'in dic'
    def get_descent(self, efoname):
         if efoname not in self.dic_all:
             print("efoname not in descent")
             return []
         return self.dic_all[efoname].return_desc_rec()
    def get_old_term(self, list_term):
        listnewterm=[]
        for term in list_term:
            if term in self.obscolescence:
               listnewterm.append(self.obscolescence[term])
        return listnewterm

       

owlobj=owl('/home/jeantristanb/Downloads/efo.owl')
writeefo=open('list.key', 'w')
writeefo.write("\n".join(owlobj.get_allkey()))
writeefo.close()
owlobj.build_tree()
listdescent=owlobj.get_descent("http://purl.obolibrary.org/obo/MONDO_0045024")
listdescent+=owlobj.get_old_term(listdescent)


writecancer=open('MONDO_0045024.descent', 'w')
writecancer.write("\n".join((list(set(listdescent)))))
writecancer.close()
