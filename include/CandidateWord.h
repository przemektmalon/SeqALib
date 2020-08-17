struct candidateWord
{
    int score;     //Score of the word
    int indexSeq1; //Index of the word on sequence 1
    int indexSeq2; //Index of the word on sequence 2
    int wordSize; //Word Size
    int maxScore; //Maximum score achieved by the alignment
    bool saved = false;

    //Overloaded equivalent operator
    candidateWord *operator=(const candidateWord b)
    {
        this->score = b.score;
        this->indexSeq1 = b.indexSeq1;
        this->indexSeq2 = b.indexSeq2;
        this->wordSize = b.wordSize;
        this->maxScore = b.maxScore;
        this->saved = b.saved;
        return this;
    }

    //TODO: Figure out which of these is actually used and provide an explanation in code itself
    //Overloaded less than operator
    bool operator<(const candidateWord &b) const
    {
        return indexSeq1 < b.indexSeq1;
    }

    //Overloaded less than or equal to operator
    bool operator<=(const candidateWord &b) const
    {
        return indexSeq2 <= b.indexSeq2;
    }
};