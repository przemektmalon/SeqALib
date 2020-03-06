struct candidateWord
{
    int score;     //Score of the word
    int indexSeq1; //Index of the word on sequence 1
    int indexSeq2; //Index of the word on sequence 2
    int wordSize;
    int maxScore; //Maximum score achieved by the alignment
    bool saved = false;

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

    bool operator<(const candidateWord &b) const
    {
        return indexSeq1 < b.indexSeq1;
    }

    bool operator<=(const candidateWord &b) const
    {
        return indexSeq2 <= b.indexSeq2;
    }
};