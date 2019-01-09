//
// Created by Mario Marchand on 17-11-28.
//


#include "ActiveSequence.h"

namespace bp = boost::python;

ActiveSequence::ActiveSequence(size_t p_nbOfSpectra, double p_window_size, bool p_vlm)
    : window_size(p_window_size), nbOfSpectra(p_nbOfSpectra), mz_avg(0), mz_lb(-50), for_vlm(p_vlm)
{
    spectraPresent.resize(nbOfSpectra);
    for (size_t i=0; i<nbOfSpectra; ++i)
    {
        spectraPresent[i] = false;
    }
}

bool ActiveSequence::empty() const
{
    return theList.empty();
}

double ActiveSequence::getAverageMz() const
{
    return mz_avg;
}

//precondition: peaks in theList must all be from different spectra
//precondition: mz_lb should have the correct value
bool ActiveSequence::isValid(const Heap & heap) const
{
    if (for_vlm)
    {
        if (theList.size() != nbOfSpectra) return false;
    }

    if(this->empty()) return false;
    if(!heap.empty())
    {
        if( std::get<2>(heap.top()) <= mz_avg*(1.0 + window_size) ) return false;
    }
    if(std::get<2>(theList.back()) > mz_avg*(1.0 + window_size)) return false;
    if(std::get<2>(theList.front()) < mz_avg*(1.0 - window_size)) return false;
    return ( mz_lb <  mz_avg*(1.0 - window_size) );
}

void ActiveSequence::advanceLowerBound()
{
    if(this->empty()) throw std::logic_error("ActiveSequence::advanceLowerBound(): ActiveSequence must be non empty");

    size_t old_size = theList.size();
    auto t = static_cast< std::tuple<size_t, size_t, double> >( theList.front() );
    theList.pop_front();

    if (std::get<0>(t) >= nbOfSpectra) throw std::logic_error("ActiveSequence::advanceLowerBound(): get<0>(t] >= nbOfSpectra");
    if(spectraPresent[std::get<0>(t)] == false) throw std::logic_error("ActiveSequence::advanceLowerBound(): that spectrum should be present");
    spectraPresent[std::get<0>(t)] = false;

    mz_lb = std::get<2>(t);

    size_t new_size = theList.size();
    if (new_size==0)
    {
        mz_avg = 0;
    }
    else
    {
        mz_avg = ((double)old_size*mz_avg - mz_lb) / (double)new_size;
    }
}

bool ActiveSequence::insert(Heap &heap, const bp::list & peak)
{
    if(heap.empty()) return false; //cannot insert one more since no more peak exists

    if(theList.empty())
    {
        std::tuple<size_t,size_t,double> t = heap.popAndFeed(peak);
        if (std::get<0>(t) >= nbOfSpectra) throw std::logic_error("ActiveSequence::insert(): get<0>(t] >= nbOfSpectra");

        if(spectraPresent[std::get<0>(t)]) throw std::logic_error("ActiveSequence::insert(): this spectra should be absent");
        spectraPresent[std::get<0>(t)] = true;
        theList.push_back(t);
        mz_avg = std::get<2>(t);
        return true;
    }
    else
    {
        std::tuple<size_t,size_t,double> t = heap.top();
        if (std::get<0>(t) >= nbOfSpectra) throw std::logic_error("ActiveSequence::insert(): get<0>(t] >= nbOfSpectra");

        size_t spectra_indx = std::get<0>(t);
        double mz = std::get<2>(t);
        if(spectraPresent[spectra_indx]) return false; //since this new peak is from a spectra already present
        size_t old_size = theList.size();
        double new_mz_avg = (old_size*mz_avg + mz) / (double)(old_size + 1);

        if (mz <= new_mz_avg * (1 + window_size))
        {
            if (std::get<2>(theList.front()) >= new_mz_avg * (1 - window_size))
            {
                //you can add this peak while the first peak in the sequence remain in the window
                spectraPresent[spectra_indx] = true;
                mz_avg = new_mz_avg;
                theList.push_back(heap.popAndFeed(peak));
                return true;
            }
        }
        return false; //since the peak is not added into the sequence
    }

}

std::ostream & operator<<(std::ostream & flux, const ActiveSequence & as)
{
    flux << "Active sequence with mz_avg = "  << as.mz_avg << " : " << std::endl;
    for(const std::tuple<size_t, size_t, double> & elem : as.theList)
    {
        auto spectra_indx = static_cast<size_t>( std::get<0>(elem) );
        auto peak_indx = static_cast<size_t>( std::get<1>(elem) );
        auto mz_value = static_cast<double>( std::get<2>(elem) );
        flux << "(" << spectra_indx << ", " << peak_indx << ", " << mz_value << ")" << std::endl;
    }
    flux << std::endl;
    return flux;
}

