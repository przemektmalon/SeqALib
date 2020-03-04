struct MUM
{
    int indexSeq1;
    int indexSeq2;
    int length;

    MUM(int newIndexSeq1, int newIndexSeq2, int newLength)
    {
        this->indexSeq1 = newIndexSeq1;
        this->indexSeq2 = newIndexSeq2;
        this->length = newLength;
    }

    MUM()
    {
        this->indexSeq1 = 0;
        this->indexSeq2 = 0;
        this->length = 0;
    }

    bool operator<(const MUM &b) const
    {
        return indexSeq1 < b.indexSeq1;
    }

    bool operator<=(const MUM &b) const
    {
        return indexSeq2 <= b.indexSeq2;
    }
};