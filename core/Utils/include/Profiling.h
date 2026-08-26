#pragma once
// =============================================================================
// Utils/Profiling.h
// -----------------------------------------------------------------------------
// Point d'entree UNIQUE pour l'instrumentation Tracy dans tout le projet.
// Corrige les conflits de redéclaration de variables en rendant chaque zone
// unique par ligne (__LINE__).
// =============================================================================

#include <string>

#ifdef TRACY_ENABLE
    #include <tracy/Tracy.hpp>

    #define PROFILE_CAT_IMPL(a, b) a##b
    #define PROFILE_CAT(a, b) PROFILE_CAT_IMPL(a, b)

    // Utilisation des versions "Named" de Tracy pour éviter le conflit de variable ___tracy_scoped_zone
    #define PROFILE_ZONE()              ZoneNamed(PROFILE_CAT(_tracy_zone_, __LINE__), true)
    #define PROFILE_ZONE_N(name)        ZoneNamedN(PROFILE_CAT(_tracy_zone_, __LINE__), name, true)
    #define PROFILE_ZONE_C(color)       ZoneNamedC(PROFILE_CAT(_tracy_zone_, __LINE__), color, true)
    #define PROFILE_ZONE_NC(name, col)  ZoneNamedNC(PROFILE_CAT(_tracy_zone_, __LINE__), name, col, true)

    // Variante "exposee" : la variable de zone garde le nom donne en premier
    // argument, ce qui permet d'appeler zone.Name(buf, len) juste apres pour
    // suffixer dynamiquement une zone partagee par plusieurs appelants (ex:
    // SystemAssembler::assemble, utilise par les 4 modules physiques) avec le
    // nom du module. Sans ca, toutes les instances templatees/partagees
    // portent le meme nom de zone statique et deviennent indistinguables
    // dans la timeline Tracy.
    #define PROFILE_ZONE_NC_VAR(varname, name, col) ZoneNamedNC(varname, name, col, true)

    // Delimite les "frames" logiques (ici : un pas de temps de simulation)
    #define PROFILE_FRAME_MARK()        FrameMark
    #define PROFILE_FRAME_MARK_N(name)  FrameMarkNamed(name)

    // Courbes numeriques
    #define PROFILE_PLOT(name, value)   TracyPlot(name, value)
    #define PROFILE_PLOT_CONFIG_INT(name) TracyPlotConfig(name, tracy::PlotFormatType::Number, true, true, 0)

    #define PROFILE_MESSAGE(txt)        TracyMessageL(txt)
#else
    #define PROFILE_ZONE()
    #define PROFILE_ZONE_N(name)
    #define PROFILE_ZONE_C(color)
    #define PROFILE_ZONE_NC(name, col)
    // CORRECTIF (bug preexistant, sans rapport avec la reproduction du
    // papier) : cette macro manquait dans la branche Tracy-desactive, ce qui
    // cassait la compilation de SystemAssembler.h (utilisee par les 4
    // modules physiques) des que TRACY_ENABLE n'est pas defini.
    #define PROFILE_ZONE_NC_VAR(varname, name, col)
    #define PROFILE_FRAME_MARK()
    #define PROFILE_FRAME_MARK_N(name)
    #define PROFILE_PLOT(name, value)
    #define PROFILE_PLOT_CONFIG_INT(name)
    #define PROFILE_MESSAGE(txt)
#endif

// -----------------------------------------------------------------------------
// Couleurs de zones (format 0xRRGGBB)[cite: 2]
// -----------------------------------------------------------------------------
#define PROFILE_COLOR_SEQUENTIAL 0xE07B39
#define PROFILE_COLOR_PARALLEL   0x3CB371
#define PROFILE_COLOR_SOLVE      0x4A90D9   // couleur "solve" generique de secours

// -----------------------------------------------------------------------------
// Couleurs par module physique.
// -----------------------------------------------------------------------------
// GenericPhysicsModule::compute_step() (Physics/include/Core/GenericPhysicsModule.h)
// est le point d'entree UNIQUE de chacun des 4 modules physiques
// (Mechanics, Polarization, Fracture, Electrostatics) a chaque iteration
// Picard. La zone parente qui l'englobe est nommee dynamiquement
// "<NomModule>::compute_step" (via zone.Name()) et coloree via
// profiling_color_for_module(), pour que le panneau "Statistics" de Tracy
// regroupe assemblage + solve + post-traitement de chaque module sous sa
// propre couleur :
//   bleu = Electrostatique, rouge = Mecanique, vert = Fracture, jaune = Polarisation
#define PROFILE_COLOR_ELECTROSTATICS 0x2E6FDE   // bleu
#define PROFILE_COLOR_MECHANICS      0xD64541   // rouge
#define PROFILE_COLOR_FRACTURE       0x27AE60   // vert
#define PROFILE_COLOR_POLARIZATION   0xF1C40F   // jaune

// -----------------------------------------------------------------------------
// Nommage des threads OpenMP pour Tracy[cite: 2].
// -----------------------------------------------------------------------------
#ifdef _OPENMP
    #include <omp.h>
    #include <cstdio>

    inline void profiling_name_omp_thread_impl() {
    #ifdef TRACY_ENABLE
        static thread_local bool s_named = false;
        if (!s_named) {
            static thread_local char s_name_buf[32];
            std::snprintf(s_name_buf, sizeof(s_name_buf), "OMP_Worker_%02d", omp_get_thread_num());
            tracy::SetThreadName(s_name_buf);
            s_named = true;
        }
    #endif
    }
    #define PROFILE_OMP_THREAD_NAME() profiling_name_omp_thread_impl()
#else
    #define PROFILE_OMP_THREAD_NAME()
#endif
