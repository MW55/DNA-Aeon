//
// Created by michael on 11/05/2021.
//

#ifndef ARITHMETIC_MODULATOR_ERROR_CORRECTION_DEBUG_LOG_H
#define ARITHMETIC_MODULATOR_ERROR_CORRECTION_DEBUG_LOG_H

#ifdef _WIN32
// make Zippy work under windows...
#define _mkdir(x) _wmkdir(reinterpret_cast<const wchar_t *>(x));
#endif

#ifndef ENABLE_DEBUG_MACRO
#define DEBUG(x)
#define VERBOSE(x)
#define INFO(x)
#else
#define DEBUG(x) do { cout << "[D] " << x << std::endl; } while (0)
#define VERBOSE(x) do { cout << "[V] " << x << std::endl; } while (0)
#define INFO(x) do { cout << "[i] " << x << endl; } while (0)
#endif
#define WARN(x) do { cout << "[!] " << x << endl; } while (0)
#define SUCCESS(x) do { cout << "[+] " << x << endl; } while (0)
#define ERROR(x) do { cerr << "[!] " << x << endl; } while (0)

#endif //ARITHMETIC_MODULATOR_ERROR_CORRECTION_DEBUG_LOG_H
