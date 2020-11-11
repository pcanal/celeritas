//----------------------------------*-C++-*----------------------------------//
// Copyright 2020 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file RootImporter.test.cc
//---------------------------------------------------------------------------//

#include "io/RootImporter.hh"
#include "io/ImportTableType.hh"
#include "io/ImportProcessType.hh"
#include "io/ImportProcess.hh"
#include "io/ImportModel.hh"
#include "physics/base/ParticleMd.hh"
#include "base/Types.hh"
#include "base/Range.hh"

#include "gtest/Main.hh"
#include "gtest/Test.hh"

using celeritas::elem_id;
using celeritas::ElementDefId;
using celeritas::GdmlGeometryMap;
using celeritas::ImportMaterial;
using celeritas::ImportMaterialState;
using celeritas::ImportModel;
using celeritas::ImportParticle;
using celeritas::ImportPhysicsTable;
using celeritas::ImportPhysicsVectorType;
using celeritas::ImportProcess;
using celeritas::ImportProcessType;
using celeritas::ImportTableType;
using celeritas::ImportVolume;
using celeritas::mat_id;
using celeritas::MaterialDefId;
using celeritas::MaterialParams;
using celeritas::MaterialParamsPointers;
using celeritas::MatterState;
using celeritas::ParticleDef;
using celeritas::ParticleDefId;
using celeritas::ParticleParams;
using celeritas::PDGNumber;
using celeritas::real_type;
using celeritas::RootImporter;
using celeritas::vol_id;

//---------------------------------------------------------------------------//
// TEST HARNESS
//---------------------------------------------------------------------------//
/*!
 * The geant-exporter-data.root is created by the app/geant-exporter using the
 * four-steel-slabs.gdml example file available in app/geant-exporter/data
 */
class RootImporterTest : public celeritas::Test
{
  protected:
    void SetUp() override
    {
        root_filename_ = this->test_data_path("io", "geant-exporter-data.root");
    }

    std::string root_filename_;
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
TEST_F(RootImporterTest, import_particles)
{
    RootImporter import(root_filename_.c_str());
    auto         data = import();

    EXPECT_EQ(19, data.particle_params->size());

    EXPECT_GE(data.particle_params->find(PDGNumber(11)).get(), 0);
    ParticleDefId electron_id = data.particle_params->find(PDGNumber(11));
    ParticleDef   electron    = data.particle_params->get(electron_id);

    EXPECT_SOFT_EQ(0.510998910, electron.mass.value());
    EXPECT_EQ(-1, electron.charge.value());
    EXPECT_EQ(0, electron.decay_constant);

    std::vector<std::string> loaded_names;
    std::vector<int>         loaded_pdgs;
    for (const auto& md : data.particle_params->md())
    {
        loaded_names.push_back(md.name);
        loaded_pdgs.push_back(md.pdg_code.get());
    }

    // clang-format off
    const std::string expected_loaded_names[] = {"gamma", "e+", "e-", "mu+",
        "mu-", "pi-", "pi+", "kaon-", "kaon+", "anti_proton", "proton",
        "anti_deuteron", "deuteron", "anti_He3", "He3", "anti_triton",
        "triton", "anti_alpha", "alpha"};
    const int expected_loaded_pdgs[] = {22, -11, 11, -13, 13, -211, 211, -321,
        321, -2212, 2212, -1000010020, 1000010020, -1000020030, 1000020030,
        -1000010030, 1000010030, -1000020040, 1000020040};
    // clang-format on

    EXPECT_VEC_EQ(expected_loaded_names, loaded_names);
    EXPECT_VEC_EQ(expected_loaded_pdgs, loaded_pdgs);
}

//---------------------------------------------------------------------------//
TEST_F(RootImporterTest, import_tables)
{
    RootImporter import(root_filename_.c_str());
    auto         data = import();

    EXPECT_GE(data.physics_tables->size(), 0);

    // Test table search
    bool lambda_kn_gamma_table = false;
    for (auto table : *data.physics_tables)
    {
        EXPECT_GE(table.physics_vectors.size(), 0);

        if (table.particle == PDGNumber{celeritas::pdg::gamma()}
            && table.table_type == ImportTableType::lambda
            && table.process == ImportProcess::compton
            && table.model == ImportModel::klein_nishina)
        {
            lambda_kn_gamma_table = true;
            break;
        }
    }
    EXPECT_TRUE(lambda_kn_gamma_table);
}

//---------------------------------------------------------------------------//
TEST_F(RootImporterTest, import_geometry)
{
    RootImporter import(root_filename_.c_str());
    auto         data = import();

    auto map = data.geometry->volid_to_matid_map();
    EXPECT_EQ(map.size(), 5);

    // Fetch a given ImportVolume provided a vol_id
    vol_id       volid  = 0;
    ImportVolume volume = data.geometry->get_volume(volid);
    EXPECT_EQ(volume.name, "box");

    // Fetch respective mat_id and ImportMaterial from the given vol_id
    mat_id         matid    = data.geometry->get_matid(volid);
    ImportMaterial material = data.geometry->get_material(matid);

    // Test material
    EXPECT_EQ(matid, 1);
    EXPECT_EQ(material.name, "G4_STAINLESS-STEEL");
    EXPECT_EQ(material.state, ImportMaterialState::solid);
    EXPECT_SOFT_EQ(material.temperature, 293.15); // [K]
    EXPECT_SOFT_EQ(material.density, 8);          // [g/cm^3]
    EXPECT_SOFT_EQ(material.electron_density,
                   2.2444324067595881e+24); // [1/cm^3]
    EXPECT_SOFT_EQ(material.atomic_density,
                   8.6993504137968536e+22);                        // [1/cm^3]
    EXPECT_SOFT_EQ(material.radiation_length, 1.7380670928095856); // [cm]
    EXPECT_SOFT_EQ(material.nuclear_int_length, 16.678055775064472); // [cm]
    EXPECT_EQ(material.elements_fractions.size(), 3);

    // Test elements within material
    const int   array_size                = 3;
    std::string elements_name[array_size] = {"Fe", "Cr", "Ni"};
    int         atomic_number[array_size] = {26, 24, 28};
    real_type   fraction[array_size]
        = {0.74621287462152097, 0.16900104431152499, 0.0847860810669534};
    real_type atomic_mass[array_size]
        = {55.845110798, 51.996130136999994, 58.693325100900005}; // [AMU]

    int i = 0;
    for (auto const& iter : material.elements_fractions)
    {
        auto elid    = iter.first;
        auto element = data.geometry->get_element(elid);

        EXPECT_EQ(element.name, elements_name[i]);
        EXPECT_EQ(element.atomic_number, atomic_number[i]);
        EXPECT_SOFT_EQ(element.atomic_mass, atomic_mass[i]);
        EXPECT_SOFT_EQ(iter.second, fraction[i]);
        i++;
    }
}

//---------------------------------------------------------------------------//
TEST_F(RootImporterTest, import_material_params)
{
    RootImporter import(root_filename_.c_str());
    auto         data = import();

    // Test material labels
    std::string material_label;
    material_label = data.material_params->id_to_label(MaterialDefId{0});
    EXPECT_EQ(material_label, "G4_Galactic");
    material_label = data.material_params->id_to_label(MaterialDefId{1});
    EXPECT_EQ(material_label, "G4_STAINLESS-STEEL");

    auto mat_host_ptr = data.material_params->host_pointers();

    // Material
    // Values differ from Geant since some are calculated at construction
    // Some differences reach ~10^-3...
    auto material = mat_host_ptr.materials[1];
    EXPECT_EQ(material.matter_state, MatterState::solid);
    EXPECT_SOFT_EQ(material.temperature, 293.15);         // [K]
    EXPECT_SOFT_EQ(material.density, 8.0080860760738588); // [g/cm^3]
    EXPECT_SOFT_EQ(material.electron_density,
                   2.2471787980801913e+24); // [1/cm^3]
    EXPECT_SOFT_EQ(material.number_density,
                   8.6993504137968536e+22);                  // [1/cm^3]
    EXPECT_SOFT_EQ(material.rad_length, 1.7363123389336561); // [cm]
    EXPECT_EQ(material.elements.size(), 3);

    // Elements within material
    // Fractions are normalized and thus may differs from the original ones
    const int array_size = 3;
    // Fe, Cr, Ni
    ElementDefId element_def_id[array_size]
        = {ElementDefId{0}, ElementDefId{1}, ElementDefId{2}};
    real_type fraction[array_size]
        = {0.74621287462152097, 0.16900104431152499, 0.0847860810669534};

    for (auto i : celeritas::range(material.elements.size()))
    {
        EXPECT_EQ(material.elements[i].element, element_def_id[i]);
        EXPECT_SOFT_EQ(material.elements[i].fraction, fraction[i]);
    }
}
