//
// Created by Mario Marchand on 17-11-28.
//

#ifndef MASSSPECTRAALIGN_MSHEAP_H
#define MASSSPECTRAALIGN_MSHEAP_H

#include <vector>
#include<tuple>
#include <algorithm>
#include <memory>
#include <iostream>

#include <boost/python.hpp>

class Heap
{

public:

    explicit Heap(const boost::python::list &);
    bool empty() const;
    size_t size() const;
    const std::tuple<size_t,size_t,double> & top() const;
    std::tuple<size_t,size_t,double> popAndFeed(const boost::python::list &);

private:

    std::vector<std::tuple<size_t, size_t, double> > theVector;

    class Comp
    {
    public:
        bool operator()(std::tuple<size_t, size_t, double> a, std::tuple<size_t, size_t, double> b) const
        {
            //implements operator< for a min_heap
            return std::get<2>(a) >= std::get<2>(b);
        }
    };

    Comp comp;

};


#endif //MASSSPECTRAALIGN_MSHEAP_H
