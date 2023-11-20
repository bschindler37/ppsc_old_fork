#ifndef _PPSC_SERIALIZATION_HPP
#define _PPSC_SERIALIZATION_HPP

// -----------------------------------------------------------------------
//
// Pseudo particle strong coupling expansion serialization utilities
//
// Author: Hugo U. R. Strand, hugo.strand@gmail.com (2016)
//
// -----------------------------------------------------------------------

#include <cereal/archives/json.hpp>
#include <cereal/types/map.hpp>

// -----------------------------------------------------------------------
namespace ppsc {
// -----------------------------------------------------------------------

typedef cereal::JSONInputArchive input_archive_type;
typedef cereal::JSONOutputArchive output_archive_type;

// -----------------------------------------------------------------------
} // end namespace ppsc
// -----------------------------------------------------------------------

#endif // _PPSC_SERIALIZATION_HPP
