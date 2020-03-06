template <typename ContainerType, typename Ty = typename ContainerType::value_type, Ty Blank = Ty(0), typename MatchFnTy = std::function<bool(Ty, Ty)>>
class MummerSA : public SequenceAligner<ContainerType, Ty, Blank, MatchFnTy>
{
  private:
    using BaseType = SequenceAligner<ContainerType, Ty, Blank, MatchFnTy>;

    ScoringSystem &Scoring = BaseType::getScoring();
    const ScoreSystemType Match = Scoring.getMatchProfit();
    const bool AllowMismatch = Scoring.getAllowMismatch();
    const ScoreSystemType Mismatch = AllowMismatch
                                         ? Scoring.getMismatchPenalty()
                                         : std::numeric_limits<ScoreSystemType>::min();

    //Build the result
    void buildResult(ContainerType &Seq1, ContainerType &Seq2, AlignedSequence<Ty, Blank> &Result)
    {

        auto &Data = Result.Data;

        //<----- Get the MUMs from a suffix tree ----->
        SuffixTree<std::string, char, '-'> suffixTree;

        std::vector<MUM> mums = suffixTree.getMUMs(Seq1, Seq2, BaseType::getMatchOperation());
        suffixTree.deleteTree();

        //<----- Sort the MUMs and perform LIS ----->

        //Sort the MUMs by index of first sequence
        std::sort(mums.begin(), mums.end());

        //Perform LIS on MUMs via index of second sequence
        if (mums.size() != 1)
        {
            //Denoting X as completeWords
            int P[mums.size()];
            int M[mums.size() + 1];

            int l = 0;
            int lo = 0;
            int hi = 0;
            for (int i = 0; i < mums.size(); i++)
            {
                //Binary search for the largest positive j <= l
                //Such that X[M[j]] <= X[i]
                lo = 1;
                hi = l;
                while (lo <= hi)
                {
                    int mid = ceil((lo + hi) / 2);
                    if (mums[M[mid]] <= mums[i])
                    {
                        lo = mid + 1;
                    }
                    else
                    {
                        hi = mid - 1;
                    }
                }

                //After searching, lo is 1 greater than the length
                //of the longest prefix of X[i]
                int newL = lo;

                //The predecessor of X[i] is the last index of
                //the subsequence of length newL - 1
                P[i] = M[newL - 1];
                M[newL] = i;

                //If we found a subsequence longer than any we've found yet
                if (newL > l)
                {
                    l = newL;
                }
            }

            //Reconstruct the longest increasing subsequence
            std::vector<MUM> newList;
            newList.resize(l);
            int k = M[l];
            for (int i = l - 1; i >= 0; i--)
            {
                newList[i] = mums[k];
                k = P[k];
            }

            mums.clear();
            mums = newList;
        }

        //<----- Bridge the gap between the MUMs----->

        //Align up to first mum
        NeedlemanWunschSA<ArrayView<ContainerType>, Ty, Blank, MatchFnTy>
            firstBridge(Scoring, BaseType::getMatchOperation());
        ArrayView<ContainerType> seq1FirstBridgeSub(Seq1);
        ArrayView<ContainerType> seq2FirstBridgeSub(Seq2);
        seq1FirstBridgeSub.sliceWindow(0, mums[0].indexSeq1);
        seq2FirstBridgeSub.sliceWindow(0, mums[0].indexSeq2);

        AlignedSequence<Ty, Blank> firstBridgeAlignment = firstBridge.getAlignment(seq1FirstBridgeSub, seq2FirstBridgeSub);
        Data.insert(Data.end(), firstBridgeAlignment.begin(), firstBridgeAlignment.end());

        //For each mum
        int i = 0;
        bool breakOut = false;
        if (mums.size() > 1)
        {
            for (i = 0; i < mums.size(); i++)
            {
                //Align mum
                //std::cout << "mumming" << std::endl;
                alignMum(mums[i], Seq1, Seq2, Result);

                //Align gap between mum ignoring overlapping mums
                int j = i + 1;
                while (true)
                {
                    if (j < mums.size())
                    {
                        //Overlapping mum so look to bridge to next mum if possible
                        if (mums[i].indexSeq1 + mums[i].length >= mums[j].indexSeq1 || mums[i].indexSeq2 + mums[i].length >= mums[j].indexSeq2)
                        {
                            //std::cout << "Overlapping" << std::endl;
                            j++;
                            continue;
                        }
                        //std::cout << "Gapping" << std::endl;
                        NeedlemanWunschSA<ArrayView<ContainerType>, Ty, Blank, MatchFnTy> bridge(Scoring, BaseType::getMatchOperation());
                        ArrayView<ContainerType> seq1Sub(Seq1);
                        ArrayView<ContainerType> seq2Sub(Seq2);
                        seq1Sub.sliceWindow(mums[i].indexSeq1 + mums[i].length, mums[j].indexSeq1);
                        seq2Sub.sliceWindow(mums[i].indexSeq2 + mums[i].length, mums[j].indexSeq2);

                        AlignedSequence<Ty, Blank> NWAlignment = bridge.getAlignment(seq1Sub, seq2Sub);
                        Data.insert(Data.end(), NWAlignment.begin(), NWAlignment.end());
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
                    //std::cout << "breaking out" << std::endl;
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

        //Align to end
        NeedlemanWunschSA<ArrayView<ContainerType>, Ty, Blank, MatchFnTy> finalBridge(Scoring, BaseType::getMatchOperation());
        ArrayView<ContainerType> seq1FinalBridgeSub(Seq1);
        ArrayView<ContainerType> seq2FinalBridgeSub(Seq2);
        seq1FinalBridgeSub.sliceWindow(mums[i - 1].indexSeq1 + mums[i - 1].length, Seq1.size());
        seq2FinalBridgeSub.sliceWindow(mums[i - 1].indexSeq2 + mums[i - 1].length, Seq2.size());

        AlignedSequence<Ty, Blank> finalBridgeAlignment = finalBridge.getAlignment(seq1FinalBridgeSub, seq2FinalBridgeSub);
        Data.insert(Data.end(), finalBridgeAlignment.begin(), finalBridgeAlignment.end());

    }

    void alignMum(MUM mum, ContainerType &Seq1, ContainerType &Seq2, AlignedSequence<Ty, Blank> &Result)
    {
        auto &Data = Result.Data;

        for (int i = 0; i < mum.length; i++)
        {
            Data.push_back(
                typename BaseType::EntryType(Seq1[i + mum.indexSeq1], Seq2[i + mum.indexSeq2], true));
        }
    }

  public:
    static ScoringSystem getDefaultScoring()
    {
        return ScoringSystem(-1, 2, -1);
    }

    MummerSA() : BaseType(getDefaultScoring(), nullptr) {}

    MummerSA(ScoringSystem Scoring, MatchFnTy Match = nullptr)
        : BaseType(Scoring, Match) {}

    virtual AlignedSequence<Ty, Blank> getAlignment(ContainerType &Seq1, ContainerType &Seq2)
    {
        AlignedSequence<Ty, Blank> Result;
        buildResult(Seq1, Seq2, Result);
        return Result;
    }
};