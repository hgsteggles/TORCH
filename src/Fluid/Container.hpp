/**
 * Provides the IntrusiveContainer and DualIntrusiveContainer classes.
 * @file Container.hpp
 *
 * @author Harrison Steggles
 *
 * @date 24/11/2014 - the first version.
 */

#ifndef CONTAINER_HPP_
#define CONTAINER_HPP_

#include <memory>
#include <iterator>
#include <type_traits>
#include <iostream>

template <class T>
class DualIntrusiveContainer;


/**
 * @class IntrusiveContainer
 *
 * @brief STL compliant container class that accepts nodes that have an IntrusiveContainer::intrusive_hook member.
 *
 * This container simply provides high performance iteration and a nice interface (hides pointers and enables range-based for looping).
 *
 * @version 0.8, 24/11/2014
 */
template <class T, class A = std::allocator<T> >
class IntrusiveContainer {
	friend DualIntrusiveContainer<T>;
public:
	typedef A allocator_type;
	typedef typename A::value_type value_type;
	typedef typename A::reference reference;
	typedef typename A::const_reference const_reference;
	typedef typename A::pointer pointer;
	typedef typename A::difference_type difference_type;
	typedef typename A::size_type size_type;

    class intrusive_hook {
        friend class ::IntrusiveContainer<T>;
        friend class DualIntrusiveContainer<T>;
    private:
        pointer next = nullptr;
    };

	class iterator {
	public:
		typedef typename A::difference_type difference_type;
		typedef typename A::value_type value_type;
		typedef typename A::reference reference;
		typedef typename A::pointer pointer;
		typedef std::forward_iterator_tag iterator_category; //or another tag

		iterator()
		: current_ptr(start_ptr)
		{ }

		iterator(pointer ptr)
		: current_ptr(ptr)
		{ }

		iterator(const iterator& it)
		: current_ptr(it.current_ptr)
		{ }

		iterator& operator=(const iterator& it) {
			current_ptr = it.current_ptr;
			return *this;
		}

		bool operator==(const iterator& it) const {
			return current_ptr == it.current_ptr;
		}
		bool operator!=(const iterator& it) const {
			return !operator==(it);
		}

		iterator& operator++() {
			current_ptr = current_ptr->hook.next;
			return *this;
		}
		reference operator*() const {
			return *current_ptr;
		}

		pointer operator->() const {
			return current_ptr;
		}

		pointer current_ptr;
	};
	class const_iterator {
	public:
		typedef typename A::difference_type difference_type;
		typedef typename A::value_type value_type;
		typedef typename A::reference const_reference;
		typedef typename A::pointer const_pointer;
		typedef std::forward_iterator_tag iterator_category; //or another tag

		const_iterator()
		: current_ptr(start_ptr)
		{ }
		const_iterator(const_pointer ptr)
		: current_ptr(ptr)
		{ }
		const_iterator (const const_iterator& it)
		: current_ptr(it.current_ptr)
		{ }
		const_iterator (const iterator& it)
		: current_ptr(it.current_ptr)
		{ }

		const_iterator& operator=(const const_iterator& it) {
			current_ptr = it.current_ptr;
			return *this;
		}

		bool operator==(const const_iterator& it) const {
			return current_ptr == it.current_ptr;
		}
		bool operator!=(const const_iterator& it) const {
			return !operator==(it);
		}

		const_iterator& operator++() {
			current_ptr = current_ptr->hook.next;
			return *this;
		}

		const_reference operator*() const {
			return *current_ptr;
		}

		const_pointer operator->() const {
			return current_ptr;
		}

		const_pointer current_ptr;
	};

	IntrusiveContainer()
	: start_ptr(nullptr)
	, end_ptr(nullptr)
	{ }

	IntrusiveContainer(pointer start, pointer end)
	: start_ptr(start)
	, end_ptr(end)
	{ }

	IntrusiveContainer(const IntrusiveContainer& c) = delete;

	IntrusiveContainer(IntrusiveContainer&& o) {
		start_ptr = o.start_ptr;
		end_ptr = o.end_ptr;
		o.start_ptr = nullptr;
		o.end_ptr = nullptr;
	}
	IntrusiveContainer<T>& operator=(IntrusiveContainer<T>&& o) {
		start_ptr = o.start_ptr;
		end_ptr = o.end_ptr;
		o.start_ptr = nullptr;
		o.end_ptr = nullptr;
		return *this;
	}

	~IntrusiveContainer() { };

	void dispose() {
		pointer ptr = start_ptr;
		while (ptr != end_ptr) {
			pointer oldptr = ptr;
			ptr = ptr->hook.next;
			delete oldptr;
			oldptr = nullptr;
		}
	}

	IntrusiveContainer<T>& operator=(const IntrusiveContainer<T>& o) = delete;

	bool operator==(const IntrusiveContainer<T>& o) const {
		return (start_ptr == o.start_ptr) && (end_ptr == o.end_ptr);
	}
	bool operator!=(const IntrusiveContainer<T>& o) const {
		return !operator==(o);
	}

	iterator begin() {
		return iterator(start_ptr);
	}
	const_iterator begin() const {
		return const_iterator(start_ptr);
	}
	const_iterator cbegin() const {
		return const_iterator(start_ptr);
	}
	iterator end() {
		return iterator(end_ptr);
	}
	const_iterator end() const {
		return const_iterator(end_ptr);
	}
	const_iterator cend() const {
		return const_iterator(end_ptr);
	}

	void swap(const IntrusiveContainer<T>& o) {
		std::swap(start_ptr, o.start_ptr);
		std::swap(end_ptr, o.end_ptr);
	}
	//size_type size();
	//size_type max_size();
	bool empty() {
		return (start_ptr == nullptr) || (start_ptr == end_ptr);
	}

	void clear() {
		start_ptr = nullptr;
		curr_ptr = nullptr;
		end_ptr = nullptr;
	}

	void push_back(pointer ptr) {
		if (ptr != nullptr) {
			if (start_ptr == nullptr)
				start_ptr = ptr;
			else
				curr_ptr->hook.next = ptr;
			curr_ptr = ptr;
		}
	}

	pointer getStartPointer() {
		return start_ptr;
	}


	pointer start_ptr = nullptr;
	pointer curr_ptr = nullptr;
	pointer end_ptr = nullptr;
};

/**
 * @class DualIntrusiveContainer
 *
 * @brief Connects and contains two IntrusiveContainers so that three types of iteration are possible.
 *
 * Constructing a DualIntrusiveContainer contains two separate IntrusiveContainers (copied from constructor arguments). Contains
 * another IntrusiveContainer which contains all the nodes. This class enables three types of iteration: over the first container;
 * over the second container; and over all the nodes contained by the first and second container.
 *
 * @see IntrusiveContainer
 *
 * @version 0.8, 24/11/2014
 */
template <class T>
class DualIntrusiveContainer {
public:

	DualIntrusiveContainer()
	{ }

	DualIntrusiveContainer(const IntrusiveContainer<T>& first, const IntrusiveContainer<T>& second)
	: m_first(first.start_ptr, first.start_ptr != nullptr ? second.start_ptr : nullptr)
	, m_second(second.start_ptr, second.end_ptr)
	, m_all(first.start_ptr != nullptr ? first.start_ptr : second.start_ptr, nullptr)
	{
		if (first.curr_ptr != nullptr)
			first.curr_ptr->hook.next = second.start_ptr;
	}

	DualIntrusiveContainer(DualIntrusiveContainer&& o) {
		std::swap(m_all, o.m_all);
		std::swap(m_first, o.m_first);
		std::swap(m_second, o.m_second);
		o.clear();
	}

	DualIntrusiveContainer<T>& operator=(DualIntrusiveContainer<T>&& o) {
		std::swap(m_all, o.m_all);
		std::swap(m_first, o.m_first);
		std::swap(m_second, o.m_second);
		o.clear();
		return *this;
	}

	void print() {
		std::cout << "dualintrusivecontainer: m_first: " << m_first.start_ptr << '\t' << m_first.end_ptr << std::endl;
		std::cout << "dualintrusivecontainer: m_second: " << m_second.start_ptr << '\t' << m_second.end_ptr << std::endl;
		std::cout << "dualintrusivecontainer: m_all: " << m_all.start_ptr << '\t' << m_all.end_ptr << std::endl;
	}

	void clear() {
		m_first.clear();
		m_second.clear();
		m_all.clear();
	}

	IntrusiveContainer<T>& getFirst() {
		return m_first;
	}

	IntrusiveContainer<T>& getSecond() {
		return m_second;
	}

	IntrusiveContainer<T>& getAll() {
		return m_all;
	}

	const IntrusiveContainer<T>& getFirst() const {
		return m_first;
	}

	const IntrusiveContainer<T>& getSecond() const {
		return m_second;
	}

	const IntrusiveContainer<T>& getAll() const {
		return m_all;
	}

	void swap(const DualIntrusiveContainer& o) {
		std::swap(m_all, o.m_all);
		std::swap(m_first, o.m_first);
		std::swap(m_second, o.m_second);
	}

private:
	DualIntrusiveContainer(const DualIntrusiveContainer&);
	DualIntrusiveContainer<T>& operator=(const DualIntrusiveContainer<T>&);

	IntrusiveContainer<T> m_first;
	IntrusiveContainer<T> m_second;
	IntrusiveContainer<T> m_all;
};


#endif /* CONTAINER_HPP_ */
