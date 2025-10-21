#include <boost/core/ignore_unused.hpp>

#ifndef ADAPT_HASHUTILS_H
#define ADAPT_HASHUTILS_H

inline void hash_combine(std::size_t &seed)
{
    boost::ignore_unused(seed);
}

template <typename T, typename... Args>
inline void hash_combine(std::size_t &seed, const T &v, Args... args)
{
    std::hash<T> hasher;
    seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    hash_combine(seed, args...);
}

template <typename T, typename... Args>
inline std::size_t hash_combine(const T &v, Args... args)
{
    boost::ignore_unused(v);
    std::size_t seed = 0;
    hash_combine(seed, args...);
    return seed;
}

template <typename Iter> std::size_t hash_range(Iter first, Iter last)
{
    std::size_t seed = 0;
    for (; first != last; ++first)
    {
        hash_combine(seed, *first);
    }
    return seed;
}

template <typename Iter>
void hash_range(std::size_t &seed, Iter first, Iter last)
{
    hash_combine(seed, hash_range(first, last));
}

struct EnumHash
{
    template <typename T> std::size_t operator()(T t) const
    {
        return static_cast<std::size_t>(t);
    }
};

//struct PairHash
//{
//    template <class T1, class T2>
//    std::size_t operator()(const std::pair<T1, T2> &p) const
//    {
//        std::size_t seed = 0;
//        auto h1          = std::hash<T1>{}(p.first);
//        auto h2          = std::hash<T2>{}(p.second);
//        hash_combine(seed, h1, h2);
//        return seed;
//    }
//};
#endif
