//
// Created by Mario Marchand on 17-11-.
//

#ifndef MASSSPECTRAALIGN_ACTIVESEQUENCE_H
#define MASSSPECTRAALIGN_ACTIVESEQUENCE_H

#include <list>
#include <vector>
#include <tuple>

#include "Heap.h"

class ActiveSequence
{
public:
    ActiveSequence(size_t nbOfSpectra, double p_window_size, bool p_vlm);
    bool isValid(const Heap &) const;
    bool empty() const;
    void advanceLowerBound();
    bool insert(Heap & heap, const boost::python::list & peak);
    double getAverageMz() const;
    friend std::ostream & operator<<(std::ostream & flux, const ActiveSequence & as);

private:
    std::list<std::tuple<size_t,size_t,double> > theList;
    std::vector<bool> spectraPresent;
    double window_size;
    size_t nbOfSpectra;
    double mz_avg;
    double mz_lb;
    bool for_vlm; //true if we want VLM identification; false if we want alignment point identification


};


#endif //MASSSPECTRAALIGN_ACTIVESEQUENCE_H
