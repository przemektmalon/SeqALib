#include <functional>
template <typename ContainerType, typename Ty = typename ContainerType::value_type, Ty Blank = Ty(0), typename MatchFnTy = std::function<bool(Ty, Ty)>>
struct Kmer
{
    size_t size;
    const ContainerType* Seq;
    size_t index;
    MatchFnTy matchFn;

    Kmer(ContainerType& sequence, size_t newSize, size_t newIndex, MatchFnTy matchFunc) : Seq(&sequence), size(newSize), index(newIndex), matchFn(matchFunc) {}

    bool operator==(const Kmer& k)
    {
        for (int i = 0; i < k.size; i++)
        {
            if (!matchFn(k.Seq[k.index + i], Seq[index + i]))
                return false;
        }

        return true;
    }

    /*bool operator>(const Kmer& k)
    {
        if (index > k.index)
            return true;
        return false;
    }*/

};