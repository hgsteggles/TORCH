#include <string>

namespace FileManagement {

std::string basename( std::string const& pathname );

void makeDirectory(const std::string& dir_name);

int makeDirectoryPath(const std::string& dir_name);

void deleteFileContents(const std::string& folder);

void copyConfigFile(const std::string& filename, const std::string& directory);

}
