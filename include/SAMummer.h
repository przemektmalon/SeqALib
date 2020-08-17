#include <iostream>
#include <fstream>
template <typename ContainerType, typename Ty = typename ContainerType::value_type, Ty Blank = Ty(0), typename MatchFnTy = std::function<bool(Ty, Ty)>>
class MummerSA : public SequenceAligner<ContainerType, Ty, Blank, MatchFnTy>
{
private:
    using BaseType = SequenceAligner<ContainerType, Ty, Blank, MatchFnTy>;

    ScoringSystem &Scoring = BaseType::getScoring();
    const ScoreSystemType Match = Scoring.getMatchProfit();
    const bool AllowMismatch = Scoring.getAllowMismatch();
    const ScoreSystemType Mismatch = AllowMismatch ? Scoring.getMismatchPenalty() : std::numeric_limits<ScoreSystemType>::min();

    // Build the result
    void buildResult(ContainerType &Seq1, ContainerType &Seq2, AlignedSequence<Ty, Blank> &Result)
    {
        auto &Data = Result.Data;

        /*if (Seq1.size() + Seq2.size() < 200)
        {
            StaticFuncs<ContainerType, Ty, Blank, MatchFnTy>::useNW(Seq1, Seq2, Result, Scoring, BaseType::getMatchOperation());
            return;
        }*/

        //<----- Get the MUMs from a suffix tree ----->
        SuffixTree<ContainerType> suffixTree;
        std::vector<MUM> mums = suffixTree.getMUMs(Seq1, Seq2, BaseType::getMatchOperation());
        suffixTree.deleteTree();

        // Use NW instead
        if (mums.size() == 0)
        {
            StaticFuncs<ContainerType, Ty, Blank, MatchFnTy>::useNW(Seq1, Seq2, Result, Scoring, BaseType::getMatchOperation());
            return;
        }

        //<----- Sort the MUMs and perform LIS ----->

        // Sort the MUMs by index of first sequence
        std::sort(mums.begin(), mums.end());

        // Perform LIS on MUMs via index of second sequence
        BaseType::longestIncreasingSubsequence(mums);


        //<----- Bridge the gap between the MUMs----->

        // Align up to first mum
        StaticFuncs<ContainerType, Ty, Blank, MatchFnTy>::bridgeNW(Seq1, Seq2, Result, Scoring, 0, 0, mums[0].indexSeq1, mums[0].indexSeq2, BaseType::getMatchOperation());

        // For each mum
        int i = 0;
        bool breakOut = false;
        if (mums.size() > 1)
        {
            for (i = 0; i < mums.size(); i++)
            {
                // Align mum
                alignMum(mums[i], Seq1, Seq2, Result);

                // Align gap between mum ignoring overlapping mums
                int j = i + 1;
                while (true)
                {
                    if (j < mums.size())
                    {
                        // Overlapping mum so look to bridge to next mum if possible
                        if (mums[i].indexSeq1 + mums[i].length > mums[j].indexSeq1 || mums[i].indexSeq2 + mums[i].length > mums[j].indexSeq2)
                        {
                            j++;
                            continue;
                        }

                        // Align the gap 
                        StaticFuncs<ContainerType, Ty, Blank, MatchFnTy>::bridgeNW(Seq1, Seq2, Result, Scoring, mums[i].indexSeq1 + mums[i].length, mums[i].indexSeq2 + mums[i].length, mums[j].indexSeq1, mums[j].indexSeq2, BaseType::getMatchOperation());

                        i = j - 1;
                        break;
                    }
                    else
                    {
                        breakOut = true;
                        break;
                    }
                }

                if (breakOut)
                {
                    i++;
                    break;
                }
            }
        }
        else
        {
            alignMum(mums[0], Seq1, Seq2, Result);
            i = 1;
        }

        // Align to end
        StaticFuncs<ContainerType, Ty, Blank, MatchFnTy>::bridgeNW(Seq1, Seq2, Result, Scoring, mums[i-1].indexSeq1 + mums[i-1].length, mums[i-1].indexSeq2 + mums[i-1].length, Seq1.size(), Seq2.size(), BaseType::getMatchOperation());
    }

    void alignMum(MUM mum, ContainerType &Seq1, ContainerType &Seq2, AlignedSequence<Ty, Blank> &Result)
    {
        auto &Data = Result.Data;

        for (int i = 0; i < mum.length; i++)
        {
            Data.push_back(typename BaseType::EntryType(
                Seq1[i + mum.indexSeq1], Seq2[i + mum.indexSeq2], true)
            );
        }
    }

public:
    static ScoringSystem getDefaultScoring() { return ScoringSystem(-1, 2, -1); }

    MummerSA() : BaseType(getDefaultScoring(), nullptr) {}

    MummerSA(ScoringSystem Scoring, MatchFnTy Match = nullptr) : BaseType(Scoring, Match) {}

    virtual AlignedSequence<Ty, Blank> getAlignment(ContainerType &Seq1, ContainerType &Seq2)
    {
        AlignedSequence<Ty, Blank> Result;
        buildResult(Seq1, Seq2, Result);
        return Result;
    }
};