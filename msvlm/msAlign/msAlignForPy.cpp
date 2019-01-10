//
// Created by Mario Marchand on 17-11-28.
//

#include <iostream>
#include <boost/python.hpp>

#include "Heap.h"
#include "ActiveSequence.h"

namespace bp = boost::python;

bp::list removeOverlaps(const std::vector<double> &tentativeAlignmentPoint, double window_size)
{
    bp::list alignmentPoint;

    size_t n = tentativeAlignmentPoint.size();
    std::vector<bool> isOverlapping(n, false);

    if (n == 0) return alignmentPoint;

    for (size_t i = 0; i < (n - 1); ++i)
    {
        if ((1.0 + window_size) * tentativeAlignmentPoint.at(i) >=
            (1.0 - window_size) * tentativeAlignmentPoint.at(i + 1))
        {
            isOverlapping[i] = true;
            isOverlapping[i + 1] = true;
        }
    }
    for (size_t i = 0; i < n; ++i)
    {
        if (!isOverlapping.at(i))
        {
            alignmentPoint.append(tentativeAlignmentPoint[i]);
        }

    }
    return alignmentPoint;
}

bp::list alignmentPointDetection(const bp::list &peak, double window_size, bool for_vlm)
{
    std::vector<double> tentativeAlignmentPoint; //aligment points that may overlap; initially empty

    auto nbOfSpectra = static_cast<size_t>(len(peak));

    Heap heap(peak); //each spectra has its first peak in the min heap h
    ActiveSequence as(nbOfSpectra, window_size, for_vlm); //initially empty (contains zero peaks)

    bool found = false;

    while (!heap.empty())
    {
        if (as.isValid(heap)) found = true;

        if (!as.insert(heap, peak)) //was not able to insert next peak in as
        {
            if (found)
            {
                tentativeAlignmentPoint.push_back(as.getAverageMz());
                found = false;
            }
            as.advanceLowerBound();
        } else //the insertion in as has succeeded
        {
            if (heap.empty()) //if there are no more peak to insert you are done
            {
                while (!as.empty())
                {
                    if (as.isValid(heap))
                    {
                        tentativeAlignmentPoint.push_back(as.getAverageMz());
                        break;
                    }
                    as.advanceLowerBound();
                }
            }
        }
    }

    return removeOverlaps(tentativeAlignmentPoint, window_size);

}


bp::list find_vlm(const bp::list  & pylistOfpylists, const bp::object & pyws)
{
    double window_size = bp::extract<double>(pyws);
    return alignmentPointDetection(pylistOfpylists, window_size, true);
}

bp::list find_alpt(const bp::list  & pylistOfpylists, const bp::object & pyws)
{
    double window_size = bp::extract<double>(pyws);
    return alignmentPointDetection(pylistOfpylists, window_size, false);

//    bp::list rl;
//    std::vector<double>::iterator it;
//    for (it = alignmentPoint.begin(); it != alignmentPoint.end(); ++it)
//    {
//        rl.append(*it);
//    }
//    return rl;
}


BOOST_PYTHON_MODULE(msAlign)
{
    bp::def("find_vlm", &find_vlm);
    bp::def("find_alpt", &find_alpt);
}
