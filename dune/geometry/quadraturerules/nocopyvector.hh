#ifndef DUNE_GEOMETRY_QUADRATURERULES_NOCOPYVECTOR_HH
#define DUNE_GEOMETRY_QUADRATURERULES_NOCOPYVECTOR_HH

#include <cstddef>

namespace Dune
{

  //! A vector of non-copyable items
  /**
   * Some items do not have a copy constructor, but do have a
   * default-constructor.  Examples include \c std::mutex, \c std::once_flag,
   * and \c std::atomic_flag.  This class can store a vector of such items.
   *
   * \note This class is an implementation detail of QuadratureRules and not
   *       part of the official Dune interface.
   */
  template<class T>
  class NoCopyVector
  {
    T* data_;
    std::size_t size_;

  public:
    //! type of the mutexes
    typedef T value_type;
    //! iterator type
    typedef T* iterator;

    //! construct
    /**
     * \param size Number of items held.
     */
    NoCopyVector(std::size_t size = 0) :
      data_(0), size_(0)
    {
      resize(size);
    }
    //! destroy
    /**
     * This will call clear()
     */
    ~NoCopyVector()
    {
      clear();
    }

    //! free the storage for all items
    /**
     * Effects: Calls destructor for all items.  Releases memory.  Invalidates
     * any iterators.
     */
    void clear() {
      if(data_)
        delete[] data_;
      data_ = 0;
      size_ = 0;
    }

    //! resize the vector
    /**
     * Effects: Calls clear().  Allocates memory for the new items.
     * Constructs the new items via their default constructor.
     *
     * \note The values of any items held by the vector prior to the call to
     *       resize() is not preserved.  After the call all items in the
     *       vector will be in the default state.
     *
     * Invalidates any iterators.
     */
    void resize(std::size_t newsize) {
      clear();
      if(newsize) {
        data_ = new T[newsize];
        size_ = newsize;
      }
    }

    //! return a reference to the i'th item
    T &operator[](std::size_t i)
    {
      return data_[i];
    }

    //! number of items currently in the vector
    std::size_t size() const
    {
      return size_;
    }

    //! iterator to the first item
    iterator begin()
    {
      return data_;
    }
    //! past-the-end iterator
    iterator end()
    {
      return data_ + size_;
    }

  };

} // namespace Dune

#endif // DUNE_GEOMETRY_QUADRATURERULES_NOCOPYVECTOR_HH
