
#include "defs.h"
#include "samples.h"
#include "logger.h"
#include <iostream>
#include <cstring>
#include <fstream>
#include <map>
#include <sstream>
#include <algorithm>
#include <limits>
// *****************************************************************************************************************
int Samples::loadSamples(vector<string> &v_samples)
{
    auto logger = LogManager::Instance().Logger();
    ok_ = true;
    uint32_t c = 0;
    for (size_t i = 0; i < v_samples.size(); i++)
    {
        if (!whichIndMap.insert(std::pair<std::string, uint32_t>(v_samples[i], c++)).second)
        {
            logger->error("Error! Duplicate sample name: {}", v_samples[i]);
            return 2;
        }
        all_samples += v_samples[i] + "\t";
    }
    all_samples = all_samples.substr(0, all_samples.length() - 1);

    no_samples = v_samples.size();
    return 0;
}
// *****************************************************************************************************************
void Samples::get_all_samples(string &str)
{
    str += all_samples;
}
// *****************************************************************************************************************
uint32_t Samples::getWhich(std::string nm)
{
    auto logger = LogManager::Instance().Logger();
    auto it = whichIndMap.find(nm);
    if (it != whichIndMap.end())
        return it->second;
    else
    {
        logger->error("There is no sample {} in the set", nm);
        ok_ = false;
        return std::numeric_limits<uint32_t>::max();
    }
}
// *****************************************************************************************************************
uint32_t *Samples::setSamples(const std::string &samples, string &str)
{
    auto logger = LogManager::Instance().Logger();
    uint32_t *smplIDs = nullptr;
    long size = 0;

    ok_ = true;
    no_samples = 0;
    if (samples[0] == '@')
    {
        std::ifstream in_samples(samples.substr(1));
        if (!in_samples.is_open())
        {
            logger->error("Error. Cannot open {} file with samples.", samples.substr(1));
            ok_ = false;
            return nullptr;
        }

        uint32_t i = 0;
        std::string item;

        while (in_samples >> item)
        {
            size++;
        }
        logger->info("Sample file size: {}", size);
        in_samples.clear();
        in_samples.seekg(0);

        smplIDs = new uint32_t[size];
        while (in_samples >> item)
        {
            str += (item + "\t");
            smplIDs[i++] = getWhich(item);
            if (!ok_)
                break;
            no_samples++;
        }
    }
    else
    {
        size = std::count(samples.begin(), samples.end(), ',') + 1;

        smplIDs = new uint32_t[size];
        char delim = ',';
        uint32_t i = 0;
        std::stringstream ss(samples);
        std::string item;
        while (getline(ss, item, delim))
        {
            str += (item + "\t");
            smplIDs[i++] = getWhich(item);
            if (!ok_)
                break;
            no_samples++;
        }
    }
    if (!ok_)
    {
        delete[] smplIDs;
        smplIDs = nullptr;
        str.clear();
        no_samples = 0;
        return nullptr;
    }

    if (!str.empty())
        str.pop_back();

    return smplIDs;
}
