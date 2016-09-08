#include "FileManagement.hpp"

#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <errno.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif /* HAVE_UNISTD_H */
#include <cstring>
#include <stdexcept>

#include <iostream>
#include <fstream>

typedef struct stat Stat;

#ifndef lint
/* Prevent over-aggressive optimizers from eliminating ID string */
const char jlss_id_mkpath_c[] = "@(#)$Id: mkpath.c,v 1.13 2012/07/15 00:40:37 jleffler Exp $";
#endif /* lint */

std::string FileManagement::basename( std::string const& pathname ) {
    return pathname.substr( pathname.find_last_of("\\/") + 1 );
}

void FileManagement::makeDirectory(const std::string& dir_name) {
	struct stat st = {0};
	if (stat(dir_name.c_str(), &st) == -1)
		mkdir(dir_name.c_str(), 0755);
}

static int do_mkdir(const char *path, mode_t mode)
{
	Stat st;
	int status = 0;

	if (stat(path, &st) != 0) {
		/* Directory does not exist. EEXIST for race condition */
		if (mkdir(path, mode) != 0 && errno != EEXIST)
			status = -1;
	}
	else if (!S_ISDIR(st.st_mode)) {
		errno = ENOTDIR;
		status = -1;
	}

	return(status);
}

int FileManagement::makeDirectoryPath(const std::string& path) {
	char *pp;
	char *sp;
	int status;
	char *copypath = strdup(path.c_str());

	status = 0;
	pp = copypath;
	while (status == 0 && (sp = strchr(pp, '/')) != 0) {
		if (sp != pp) {
			/* Neither root nor double slash in path */
			*sp = '\0';
			status = do_mkdir(copypath, 0755);
			*sp = '/';
		}
		pp = sp + 1;
	}
	if (status == 0)
		status = do_mkdir(path.c_str(), 0755);
	free(copypath);

	return (status);
}

void FileManagement::deleteFileContents(const std::string& folder){
	struct dirent *next_file;
	DIR *dir; // These are data types defined in the "dirent" header.
	char filepath[256];
	dir = opendir(folder.c_str() );
	if (!dir) {
		throw std::runtime_error("deleteFileContents: directory " + folder + " does not exist.");
	}
	else {
		while ((next_file = readdir(dir) )) {
			// Build the full path for each file in the folder.
			sprintf(filepath, "%s/%s", folder.c_str(), next_file->d_name);
			remove(filepath);
		}
	}

	closedir(dir);
}

void FileManagement::copyConfigFile(const std::string& filename, const std::string& directory) {
	std::ifstream  src(filename, std::ios::binary);
	std::ofstream  dst(directory + "/" + basename(filename), std::ios::binary);
	dst << src.rdbuf();
}
