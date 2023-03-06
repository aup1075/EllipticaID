#include "nsns_header.h"
#include "id_reader_lib.h"
#include "elliptica_id_reader_lib.h"

/* prefix parameters came from evo codes, should be lower case  */
#define EVO_ "evo_"

#define STR_LEN_MAX (9999)

void nsns_export_id_generic(Elliptica_ID_Reader_T *const idr);
void nsns_set_evo_fields_generic(Grid_T *const grid);

