/*
SUPER CONFIDENTIAL COPYRIGHT TRADEMARK PATENT BY PREM MALON
ANY VIOLATIONS WILL BRING THE DEATH PENALTY
*/

#include <vector>

template<typename T>
class AlignVec {
public:
    AlignVec()
    {
        AlignVec(32);
    }
    ~AlignVec() = default;
    AlignVec(const AlignVec& AV)
    {
        m_container = AV.m_container;
        left = AV.left;
        right = AV.right;
    }

    AlignVec(AlignVec&&) = delete;
    AlignVec& operator=(AlignVec&&) = delete;
    AlignVec& operator=(const AlignVec& AV)
    {
        m_container = AV.m_container;
        left = AV.left;
        right = AV.right;
    }

    AlignVec(size_t initialSize) {
        m_container.resize(initialSize);
        left = initialSize / 2;
        right = left + 1;
        assert(right < initialSize && "Bad initial size");
        assert(left >= 0 && "Bad initial size");
    }

    void push_front(const T& val) {
        m_container[left] = val;
        if (--left < 0)
            grow();
    }

    void push_back(const T& val) {
        m_container[right] = val;
        if (++right >= m_container.size())
            grow();
    }

    void push_back(const AlignVec& V)
    {
        for (int i = 0; i < V.size(); i++)
        {
            push_back(V[i]);
        }

    }

    auto begin()
    {
        return m_container.begin() + left + 1;
    }

    auto end()
    {
        return m_container.end() + right - 1;
    }

    auto size()
    {
        return m_container.size();
    }

    void clear()
    {
        m_container.clear();
        AlignVec(32);
    }

    /*template<class InputIt>
    void insert(InputIt thisIter, InputIt first, InputIt end)
    {

    }*/

    T& operator[](size_t index) {
        assert(1 + left + index < right);
        return m_container[1 + left + index];
    }

private:
    std::vector<T> m_container;

    size_t left, right;  // left occupies the left of the leftmost populated index
                         // right occupies the right of the rightmost populated index

    void grow() {
        // CRUDE WAY OF DOING THIS! You can figure out the better, in place, way ;)
        std::vector<T> temp{ m_container };
        m_container.resize(m_container.size() * 1.5);
        m_container.clear();
        size_t mid = m_container.size() / 2;
        left = mid - (temp.size() / 2);
        right = left + temp.size();

        // THIS IS ALSO CRUDE! But I'm working rn so. Theres a better way.
        while (right - left < temp.size() + 2) {
            ++right;
        }
        while (right - left > temp.size() + 2) {
            --right;
        }
        std::memcpy(static_cast<void*>(&m_container[left + 1]), static_cast<void*>(temp.data()), sizeof(T) * temp.size());
    }

};
