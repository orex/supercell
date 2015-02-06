/* 
 * File:   archive_extent_patch.h
 * Author: kirill
 *
 * Created on February 5, 2015, 2:48 PM
 */

int cmpsuff(const char *str, const char *suffix)
{
  size_t length_str     = strlen(str);
  size_t length_suffix  = strlen(suffix);
  
  return strcmp((str + length_str - length_suffix), suffix);
}

int archive_write_set_format_filter_by_ext(struct archive *a, const char *filename)
{
  static struct { const char *name; int (*format)(struct archive *); int (*filter)(struct archive *);  } names[] =
  {
	{ ".7z",	archive_write_set_format_7zip,    archive_write_add_filter_none},
	{ ".zip",	archive_write_set_format_zip,     archive_write_add_filter_none},
	{ ".jar",	archive_write_set_format_zip,     archive_write_add_filter_none},
	{ ".a",	        archive_write_set_format_ar_bsd,  archive_write_add_filter_none},
	{ ".ar",	archive_write_set_format_ar_bsd,  archive_write_add_filter_none},
	{ ".cpio",	archive_write_set_format_cpio,    archive_write_add_filter_none},
	{ ".shar",	archive_write_set_format_shar,    archive_write_add_filter_none},
	{ ".iso",	archive_write_set_format_iso9660, archive_write_add_filter_none},
	{ ".tar",	archive_write_set_format_gnutar,  archive_write_add_filter_none},
	{ ".tgz",	archive_write_set_format_gnutar,  archive_write_add_filter_gzip},
	{ ".tar.gz",	archive_write_set_format_gnutar,  archive_write_add_filter_gzip},
	{ ".tar.bz2",	archive_write_set_format_gnutar,  archive_write_add_filter_bzip2},
	{ ".tar.xz",	archive_write_set_format_gnutar,  archive_write_add_filter_xz},
	{ NULL,		NULL,                             NULL }
  };
  
  int i;

  for (i = 0; names[i].name != NULL; i++) 
  {
    if (cmpsuff(filename, names[i].name) == 0)
    { 
      if( (names[i].format)(a) == ARCHIVE_OK )
        return ((names[i].filter)(a));
      else
        return ARCHIVE_FATAL;
      }       
  }

  archive_set_error(a, EINVAL, "No such format '%s'", filename);
  //a->state = ARCHIVE_STATE_FATAL;
  return (ARCHIVE_FATAL);
}




