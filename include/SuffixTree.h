template <typename T_Key, typename T_Value>
class SuffixTree
{

  public:
    SuffixTree()
    {
    }

    //Insert into the suffix tree
    void insert(const T_Key &key, const T_Value &value)
    {
        terminatedKey = key + "$";
        const auto start = std::cbegin(terminatedKey);
        const auto end = std::next(start, terminatedKey.size());
        const InternalKey substr(start, end);

        if (str_empty(substr))
        {
            return;
        }

        insertUkkonen(substr, value);
    }

    //Find the values for any key containing this substring
    std::set<T_Value> find(const T_Key &substr)
    {
        std::set<T_Value> result;

        if (substr.empty())
        {
            return result;
        }

        const InternalKey substr_key(std::cbegin(substr), std::cend(substr));

        return find(substr_key);
    }

    class InternalKey
    {

      public:
        InternalKey(const KeyIterator &start)
    }
}