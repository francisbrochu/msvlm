//
// Created by Mario Marchand on 17-11-28.
//

#include "Heap.h"

namespace bp = boost::python;

Heap::Heap(const bp::list & peak)
{
    for(ssize_t i=0 ; i<len(peak); ++i)
    {
        bp::list spectra = bp::extract<bp::list>(peak[i]);
        if(len(spectra)<=0) throw std::logic_error("Heap::Heap(): a spectrum contains no peaks");
        double spectra_peak = bp::extract<double>(spectra[0]);
        theVector.emplace_back( std::tuple<size_t,size_t,double>(i, 0, spectra_peak) );
    }
    make_heap(theVector.begin(), theVector.end(), comp);
}

bool Heap::empty() const
{
    return theVector.empty();
}
size_t Heap::size() const
{
    return theVector.size();
}


const std::tuple<size_t, size_t, double> & Heap::top() const
{
    if(theVector.empty()) throw std::logic_error("Heap::top(): the heap must be non empty");
    return theVector[0];
}

std::tuple<size_t, size_t, double> Heap::popAndFeed(const bp::list & peak)
{
    if(theVector.empty()) throw std::logic_error("Heap::popAndFeed(): theVector must be non empty");

    auto returned_peak = static_cast< std::tuple<size_t, size_t, double> >( theVector.at(0) );

    size_t spectra_indx = std::get<0>(returned_peak);
    size_t peak_indx = std::get<1>(returned_peak) + 1; //must point to the next available peak from same spectra

    pop_heap(theVector.begin(), theVector.end(), comp); //swap 1st and last element and reconstruct heap without last element

    bp::list spectra = bp::extract<bp::list>(peak[spectra_indx]);

    if( peak_indx < static_cast<size_t>(len(spectra)) )
    {
        theVector[theVector.size() - 1] = std::tuple<size_t, size_t, double>(spectra_indx, peak_indx,
                                                                        bp::extract<double>(spectra[peak_indx]));
        push_heap(theVector.begin(), theVector.end(), comp);
    }
    else
    {
        theVector.pop_back();
    }

    return returned_peak;
}

