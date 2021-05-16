/*
 * File:   archive_extent_patch.h
 * Author: kirill
 *
 * Created on February 5, 2015, 2:48 PM
 */

static
int cmpsuff(const char *str, const char *suffix)
{
  size_t length_str, length_suffix;

  if ((str == NULL) || (suffix == NULL))
    return -1;

  length_str = strlen(str);
  length_suffix = strlen(suffix);

  if (length_str >= length_suffix) {
    return strcmp(str + (length_str - length_suffix), suffix);
  } else {
    return -1;
  }
}

static const
struct { const char *name; int (*format)(struct archive *); int (*filter)(struct archive *);  } names[] =
    {
        { ".7z",	archive_write_set_format_7zip,            archive_write_add_filter_none},
        { ".zip",	archive_write_set_format_zip,             archive_write_add_filter_none},
        { ".jar",	archive_write_set_format_zip,             archive_write_add_filter_none},
        { ".cpio",	archive_write_set_format_cpio,            archive_write_add_filter_none},
        { ".iso",	archive_write_set_format_iso9660,         archive_write_add_filter_none},
        { ".a",	        archive_write_set_format_ar_svr4,         archive_write_add_filter_none},
        { ".ar",	archive_write_set_format_ar_svr4,         archive_write_add_filter_none},
        { ".tar",	archive_write_set_format_pax_restricted,  archive_write_add_filter_none},
        { ".tgz",	archive_write_set_format_pax_restricted,  archive_write_add_filter_gzip},
        { ".tar.gz",	archive_write_set_format_pax_restricted,  archive_write_add_filter_gzip},
        { ".tar.bz2",	archive_write_set_format_pax_restricted,  archive_write_add_filter_bzip2},
        { ".tar.xz",	archive_write_set_format_pax_restricted,  archive_write_add_filter_xz},
        { NULL,		NULL,                             NULL }
    };

static int get_array_index(const char *name)
{
  int i;

  for (i = 0; names[i].name != NULL; i++)
  {
    if (cmpsuff(name, names[i].name) == 0)
      return i;
  }
  return -1;

}

int
archive_write_set_format_filter_by_ext(struct archive *a, const char *filename)
{
  int names_index = get_array_index(filename);

  if (names_index >= 0)
  {
    int format_state = (names[names_index].format)(a);
    if (format_state == ARCHIVE_OK)
      return ((names[names_index].filter)(a));
    else
      return format_state;
  }

  archive_set_error(a, EINVAL, "No such format '%s'", filename);
  //a->state = ARCHIVE_STATE_FATAL;
  return (ARCHIVE_FATAL);
}

