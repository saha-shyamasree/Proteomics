
class AminoAcidVariation:
    """This clas holds aminoacid variations from blast alignment"""
    def __init__(self,subjectId, queryId, start, refSeq, altSeq, Type, chrm, vid, info):
        self.subjectId=subjectId
        self.queryId=queryId
        self.pos=start
        self.ref=refSeq
        self.alt=altSeq
        self.type=Type
        self.chr=chrm
        self.vid=vid
        self.alingmentInfo=info
    def getSubjectId(self):
        return self.subjectId
    def getQueryId(self):
        return self.queryId
    def getPos(self):
        return self.pos
    def getREF(self):
        return self.ref
    def getALT(self):
        return self.alt
    def getType(self):
        return self.type
    def setSubjectId(self,subjectId):
        self.subjectId=subjectId
    def setQueryId(self,queryId):
        self.queryId=queryId
    def setPos(self,pos):
        self.pos=pos
    def setREF(self,ref):
        self.ref=ref
    def setALT(self,alt):
        self.alt=alt
    def setALT(self,Type):
        self.type=Type
    def toString(self):
        if int(self.alingmentInfo['subjectStart'])!=1:
            prefix=int(self.alingmentInfo['subjectStart'])-1;
            return self.chr+"\t"+str(self.pos+prefix)+"\t"+str(self.alingmentInfo['qSerialNo'])+"."+str(self.vid)+"\t"+self.ref+"\t"+self.alt+"\t"+"SubjectId="+self.subjectId+";QueryId="+self.queryId+";Alignment=[QueryLength="+self.alingmentInfo['qLen']+":QueryStart="+self.alingmentInfo['queryStart']+":QueryEnd="+self.alingmentInfo['queryEnd']+":SubjectLength="+self.alingmentInfo['sLen']+":SubjectStart="+self.alingmentInfo['subjectStart']+":SubjectEnd="+self.alingmentInfo['subjectEnd']+"];Type:"+self.type;
        else:
            return self.chr+"\t"+str(self.pos)+"\t"+str(self.alingmentInfo['qSerialNo'])+"."+str(self.vid)+"\t"+self.ref+"\t"+self.alt+"\t"+"SubjectId="+self.subjectId+";QueryId="+self.queryId+";Alignment=[QueryLength="+self.alingmentInfo['qLen']+":QueryStart="+self.alingmentInfo['queryStart']+":QueryEnd="+self.alingmentInfo['queryEnd']+":SubjectLength="+self.alingmentInfo['sLen']+":SubjectStart="+self.alingmentInfo['subjectStart']+":SubjectEnd="+self.alingmentInfo['subjectEnd']+"];Type:"+self.type;
    def printVar(self):
        print(self.chr+"\t"+str(self.pos)+"\t"+str(self.vid)+"\t"+self.ref+"\t"+self.alt+"\t"+"ProtId:"+self.subjectId+";ORFId:"+self.queryId+";Type:"+self.type)
    @staticmethod
    def printHeader():
        return "#Chr\tPOS within Protein\tID\tREF\tALT\tINFO"
    
    
    
