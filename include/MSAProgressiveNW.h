template <typename ContainerType, typename Ty = typename ContainerType::value_type, Ty Blank = Ty(0), typename MatchFnTy = std::function<bool(Ty, Ty)>>
class NeedlemanWunschMSA : public SequenceAligner<ContainerType, Ty, Blank, MatchFnTy>
{
private:

    ScoreSystemType* Matrix;
    size_t NumRows;
    size_t NumCols;

    using BaseType = SequenceAligner<ContainerType, Ty, Blank, MatchFnTy>;

    void alignSeq(ContainerType& Seq1, ContainerType& Seq2, AlignedSequence<Ty, Blank>& Result)
    {
        StaticFuncs<ContainerType, Ty, Blank, MatchFnTy>::useNW(Seq1, Seq2, Result, BaseType::getScoring(), BaseType::getMatchOperation());
    }

    // Align a Sequence with an Alignment
    void alignSeqAln(ContainerType& Seq, AlignedSequence<Ty, Blank>& Aln)
    {
        const size_t NumRows = Seq.size();
        const size_t NumCols = Aln.size();

        Matrix = new int[NumRows * NumCols];

        ScoringSystem& Scoring = BaseType::getScoring();
        const ScoreSystemType Gap = Scoring.getGapPenalty();
        const ScoreSystemType Match = Scoring.getMatchProfit();
        const bool AllowMismatch = Scoring.getAllowMismatch();
        const ScoreSystemType Mismatch = AllowMismatch ? Scoring.getMismatchPenalty() : std::numeric_limits<ScoreSystemType>::min();


        for (unsigned i = 0; i < NumRows; i++)
            Matrix[i * NumCols + 0] = i * Gap;
        for (unsigned j = 0; j < NumCols; j++)
            Matrix[0 * NumCols + j] = j * Gap;


        if (AllowMismatch)
        {
            for (unsigned i = 1; i < NumRows; i++)
            {
                for (unsigned j = 1; j < NumCols; j++)
                {
                    // Calculate substitution score using all the alignment characters
                    ScoreSystemType Similarity = BaseType::match(Seq[i-1], Aln[j-1].get(0)) ? Match : Mismatch;
                    for (unsigned k = 1; k < Aln.nObjects; k++)
                    {
                        Similarity += BaseType::match(Seq[i - 1], Aln[j - 1].get(k)) ? Match : Mismatch;
                    }

                    Similarity = Similarity / (NumRows * NumCols);

                    ScoreSystemType Diagonal = Matrix[(i - 1) * NumCols + j - 1] + Similarity;
                    ScoreSystemType Upper = Matrix[(i - 1) * NumCols + j] + Gap;
                    ScoreSystemType Left = Matrix[i * NumCols + j - 1] + Gap;
                    ScoreSystemType Score = std::max(std::max(Diagonal, Upper), Left);
                    Matrix[i * NumCols + j] = Score;

                }
            }
        }
        else
        {
            for (unsigned i = 0; i < NumRows; i++)
            {
                for (unsigned j = 0; j < NumRows; j++)
                {
                    // Either the sequence matches every character in the alignment column and is given a matching score
                    // Or not and we align with a gap
                    // However, since we know that we are not allowing mismatches then it's safe to say that
                    // The alignment does not contain any mismatches from it's previous builds

                    ScoreSystemType Diagonal = BaseType::match(Seq[i-1], Aln[j-1].get(0)) ? (Matrix[(i - 1) * NumCols + j - 1] + Match) : Mismatch;
                    ScoreSystemType Upper = Matrix[(i - 1) * NumCols + j] + Gap;
                    ScoreSystemType Left = Matrix[i * NumCols + j - 1] + Gap;
                    ScoreSystemType Score = std::max(std::max(Diagonal, Upper), Left);
                    Matrix[i * NumCols + j] = Score;
                }
            }
        }

    }

    // Build the resulting alignment
    void buildResult(ContainerType& Seq, AlignedSequence<Ty, Blank>& Aln)
    {
        AlignedSequence<Ty, Blank> Result;

        ScoringSystem& Scoring = BaseType::getScoring();
        const ScoreSystemType Gap = Scoring.getGapPenalty();
        const ScoreSystemType Match = Scoring.getMatchProfit();
        const bool AllowMismatch = Scoring.getAllowMismatch();
        const ScoreSystemType Mismatch = AllowMismatch ? Scoring.getMismatchPenalty() : std::numeric_limits<ScoreSystemType>::min();

        int i = NumRows - 1, j = NumCols - 1;

        while (i > 0 || j > 0)
        {
            if (i > 0 && j > 0)
            {
                bool matchesBottom = BaseType::match(Seq[i], Aln.Data[j].get(Aln.nObjects-1));

                ScoreSystemType Score = std::numeric_limits<ScoreSystemType>::min();

                if (AllowMismatch)
                {

                    Score = matchesBottom ? Match : Mismatch;
                    for (unsigned k = 0; k < Aln.nObjects-1; k++)
                    {
                        Score += Aln.Data[j].getMatch(k) ? Match : Mismatch;
                    }

                    Score = Score / (NumRows * NumCols);

                    Score = Matrix[(i - 1) * NumCols + j - 1] + Score;
                }
                else
                {
                    Score = matchesBottom ? (Matrix[(i - 1) * NumCols + j - 1] + Match) : Mismatch;
                }

                if (Matrix[i * NumCols + j] == Score)
                {
                    if (matchesBottom || AllowMismatch)
                    {
                        Result.Data.push_front(
                            typename BaseType::EntryType(Aln.Data[j - 1], Seq[i - 1], matchesBottom)
                        );

                    }

                    i--;
                    j--;
                    continue;
                }
            }
            if (i > 0 && Matrix[i * NumCols + j] == (Matrix[(i - 1) * NumCols + j] + Gap))
            {
                // Align Seq with gaps
                Aln.insert(Seq[i - 1], j);
                i--;
            }
            else
            {
                // Do Nothing?
                j--;
            }
        }

        Aln = Result;
    }

    // Align an Alignment with an Alignment
    void alignAlnAln(AlignedSequence<Ty, Blank>& Aln1, AlignedSequence<Ty, Blank>& Aln2)
    {

    }

public:
    static ScoringSystem getDefaultScoring()
    {
        return ScoringSystem(-1, 2, -1);
    }

    NeedlemanWunschMSA() : BaseType(getDefaultScoring(), nullptr), Matrix(nullptr){}

    NeedlemanWunschMSA(ScoringSystem Scoring, MatchFnTy Match = nullptr) : BaseType(Scoring, Match), Matrix(nullptr){}

    virtual AlignedSequence<Ty, Blank> getAlignment(ContainerType& Seq1, ContainerType& Seq2)
    {
        AlignedSequence<Ty, Blank> Result;
        StaticFuncs<ContainerType, Ty, Blank, MatchFnTy>::useNW(Seq1, Seq2, Result, BaseType::getScoring(), BaseType::getMatchOperation());
        return Result;
    }

    virtual AlignedSequence<Ty, Blank> getAlignment(std::vector<ContainerType> Seqs)
    {
        assert(Seqs.size() >= 2 && "Number of sequences entered too small!");

        if (Seqs.size() == 2)
            return getAlignment(Seqs[0], Seqs[1]);

        // DO UPGMA OR SOMETHING
        // For now align each sequence progressively to the current alignment
        // Though should test aligning aln to aln

        // Align Seq1 with Seq2
        AlignedSequence<Ty, Blank> Result;
        StaticFuncs<ContainerType, Ty, Blank, MatchFnTy>::useNW(Seqs[0], Seqs[1], Result, BaseType::getScoring(), BaseType::getMatchOperation());

        // Align each extra sequence to this first alignment
        for (int i = 2; i < Seqs.size(); i++)
        {
            alignSeqAln(Seqs[i], Result);
            buildResult(Seqs[i], Result);
        }

        return Result;
    }
};
