#pragma once

#include "Polynomial.h"
#include <experimental/filesystem>

namespace Check {
	void isDirectory(std::string dir) {
		if (!std::experimental::filesystem::is_directory(dir)) {
			cout << dir << " is not a directory";
			exit(1);
		}
	}

	void isPathToFile(std::string path_to_file) {
		if (!std::experimental::filesystem::exists(path_to_file)) {
			cout << path_to_file << " not found";
			exit(1);
		}
	}

}
