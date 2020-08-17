#include <functional>

// Exists to provide functionality to algorithms which make use of other algorithms
// Such as providing bridging with NW 
// Cannot simply be inherited from the SequenceAligner class
template <typename ContainerType, typename Ty = typename ContainerType::value_type, Ty Blank = Ty(0), typename MatchFnTy = std::function<bool(Ty, Ty)>>
class StaticFuncs
{
public:
    // When something goes wrong in an algorithm or it is determined that using NW would be more beneficial
    // Then run this which runs NW instead
    void static useNW(ContainerType& Seq1, ContainerType& Seq2, AlignedSequence<Ty, Blank>& Result, ScoringSystem Scoring, MatchFnTy match)
    {
        auto& Data = Result.Data;
        Data.clear();

        NeedlemanWunschSA<ArrayView<ContainerType>, Ty, Blank, MatchFnTy> NW(Scoring, match);
        ArrayView<ContainerType> seq1Sub(Seq1);
        ArrayView<ContainerType> seq2Sub(Seq2);
        seq1Sub.sliceWindow(0, Seq1.size());
        seq2Sub.sliceWindow(0, Seq2.size());
        AlignedSequence<Ty, Blank> NWAlignment = NW.getAlignment(seq1Sub, seq2Sub);
        Data.insert(Data.end(), NWAlignment.Data.begin(), NWAlignment.Data.end());
    }

    //Bridging using NW to align the spaces between the anchors found in algorithms like BLAST or MUMmer
    void static bridgeNW(ContainerType& Seq1, ContainerType& Seq2, AlignedSequence<Ty, Blank>& Result, ScoringSystem Scoring, int idx1, int idx2, int endIdx1, int endIdx2, MatchFnTy match)
    {
        auto& Data = Result.Data;

        NeedlemanWunschSA<ArrayView<ContainerType>, Ty, Blank, MatchFnTy> NWBridge(Scoring, match);
        ArrayView<ContainerType> seq1BridgeSub(Seq1);
        ArrayView<ContainerType> seq2BridgeSub(Seq2);
        seq1BridgeSub.sliceWindow(idx1, endIdx1);
        seq2BridgeSub.sliceWindow(idx2, endIdx2);

        AlignedSequence<Ty, Blank> bridgeAlignment = NWBridge.getAlignment(seq1BridgeSub, seq2BridgeSub);
        Data.insert(Data.end(), bridgeAlignment.Data.begin(), bridgeAlignment.Data.end());
    }
};