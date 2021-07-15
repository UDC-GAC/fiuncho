#include <cstring>
#include <memory>

template <class T> class StaticStack
{
  StaticStack(const StaticStack<T> &) = delete;
  StaticStack(StaticStack<T> &&) = default;

  public:
    StaticStack(size_t item_size, size_t max_items)
        : item_size(item_size),
          array_ptr(std::make_unique<char[]>(item_size * max_items)),
          array((char *)array_ptr.get()), top(0)
    {
    }

    inline void push(const T &item)
    {
        memcpy(array + top++ * item_size, &item, item_size);
    }

    inline void pop(T &item)
    {
        memcpy(&item, array + --top * item_size, item_size);
    }

    inline size_t size() { return top; }

    inline bool empty() { return top == 0; }

  public:
    size_t item_size;
    std::unique_ptr<char[]> array_ptr;
    char *array;
    size_t top;
};