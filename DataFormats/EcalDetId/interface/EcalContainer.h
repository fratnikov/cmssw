#ifndef ECALDETID_ECALCONTAINER_H
#define ECALDETID_ECALCONTAINER_H

#include "DataFormats/DetId/interface/DetId.h"
#include "FWCore/Utilities/interface/Exception.h"
//#include <type_traits> //C++11
#include <boost/type_traits/conditional.hpp>
#include <vector>
#include <utility>
#include <algorithm>

//#include <iostream> 



/* a generic container for ecal items
 * provides access by hashedIndex and by DetId...
 */

template<typename DetId, typename T>
class EcalContainer {

        public:

                typedef EcalContainer<DetId, T> self;
                typedef T Item;
                typedef Item value_type;
                typedef typename std::vector<Item> Items; 
                typedef typename std::vector<Item>::const_iterator const_iterator;
                typedef typename std::vector<Item>::iterator iterator;
		typedef typename boost::conditional<DetId::kSizeForDenseIndexing<((uint32_t) -1), uint32_t, uint64_t>::type Index32;
		typedef typename boost::conditional<DetId::kSizeForDenseIndexing<((uint16_t) -1), uint16_t, Index32>::type Index16;
		typedef typename boost::conditional<DetId::kSizeForDenseIndexing<((uint8_t) -1), uint8_t, Index16>::type Index;

                   
                EcalContainer() {
		  checkAndResize();
		  m_items.push_back (Item ()); // reserve index 0
		}

                void insert(std::pair<uint32_t, Item> const &a) {
                        (*this)[a.first] = a.second;
                }

                inline const Item & item(size_t hashid) const {
		  Index index = m_refs[hashid];
		  if (index == 0) return  m_items[0]; // not defined
		  return m_items[index];
                }

                inline const Items & items() const {
                        return m_items;
                }

                inline Item & operator[](uint32_t rawId) {
		  checkAndResize();
		  DetId id(rawId);
		  if ( !isValidId(id) ) {
		     throw cms::Exception("No support") << "EcalContainer-> operator[] is called with non ECAL DetId";
		  }
		  size_t hashIndex = id.hashedIndex();
		  Index index = m_refs[hashIndex];
		  if (index == 0) { // not defined, insert
		    m_items.push_back (Item());
		    index = m_refs[hashIndex] = m_items.size()-1;
		  }
		  return m_items[index];
                }


		void checkAndResize() {
		  if (m_refs.size()==0) {
		    //		    std::cout << "resizing to " << DetId::kSizeForDenseIndexing << std::endl;
		    m_refs.resize(DetId::kSizeForDenseIndexing, 0);
		  }
		}


		void checkAndResize( size_t priv_size ) {
		  // this method allows to resize the vector to a specific size forcing a specific value
		  if (m_refs.size()==0) {
		  // F.R. this should never be called as EcalContainer is initiated in the constructor
		    //		    std::cout << "resizing to " << priv_size << std::endl;
		    m_refs.resize(priv_size, 0);
		  }
		  else {
		    throw cms::Exception("No support") << "EcalContainer-> attempt to resize non empty container";
		  }
		}

                inline Item const & operator[](uint32_t rawId) const {
		  //                        if (m_items.size()==0) {
		  //	  std::cout << "resizing to " << DetId::kSizeForDenseIndexing << std::endl;
                  //              m_items.resize((size_t) DetId::kSizeForDenseIndexing);
                  //      }
                        static Item dummy;
                        DetId id(rawId);
                        if ( !isValidId(id) ) return m_items[0];
			Index index = m_refs[id.hashedIndex()];
			if (index == 0) return m_items[0];
                        return m_items[index];
                }

                inline const_iterator find(uint32_t rawId) const {
                        DetId ib(rawId);
                        if ( !isValidId(ib) ) return m_items.end();
			Index index = m_refs[ib.hashedIndex()];
			if (index == 0) return  m_items.end(); // not defined
                        return m_items.begin() + index;
                }

                inline const_iterator begin() const {
                        return m_items.begin();
                }

                inline const_iterator end() const {
                        return m_items.end();
                }

                inline size_t size() const {
                        return m_items.size();
                }

        private:

                // not protected on EB <--> EE swap -- FIXME?
                inline bool isValidId(const DetId id) const {
                        return id.det() == ::DetId::Ecal;
                };

                std::vector<Item> m_items;
		std::vector<Index> m_refs;

};



#endif // ECALCONTAINER
