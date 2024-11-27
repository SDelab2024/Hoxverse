// src/HoxGeneratorApp.js
import React, { useState, useEffect } from 'react';
import ScatterPlot from './components_custom/ScatterPlot';
import DataGenerator from './components_custom/DataGenerator';
import LabelSelector from './components_custom/LabelSelector';
import axios from 'axios';

function HoxGeneratorApp() {
  const [generatedDatasets, setGeneratedDatasets] = useState([]);
  const [selectedLabels, setSelectedLabels] = useState([]);

  useEffect(() => {
    const storedDatasets = JSON.parse(localStorage.getItem('datasets') || '[]');
    setGeneratedDatasets(storedDatasets);
  }, []);

  const handleDataGenerated = (label) => {
    const newDataset = {
      label: label,
      data: [],
      backgroundColor: `rgba(${Math.floor(Math.random() * 255)}, ${Math.floor(Math.random() * 255)}, ${Math.floor(Math.random() * 255)}, 1)`,
      csvFile: `${process.env.REACT_APP_BACKEND_URL}/data/${label}.csv`,
    };

    axios
      .get(newDataset.csvFile)
      .then((response) => {
        newDataset.data = response.data;
        const updatedDatasets = [...generatedDatasets, newDataset];
        setGeneratedDatasets(updatedDatasets);
        localStorage.setItem('datasets', JSON.stringify(updatedDatasets));
      })
      .catch((error) => {
        console.error('There was an error fetching the data!', error);
      });
  };

  const handleLabelSelected = (label) => {
    const index = selectedLabels.indexOf(label);
    if (index === -1) {
      setSelectedLabels([...selectedLabels, label]);
    } else {
      const updatedLabels = [...selectedLabels];
      updatedLabels.splice(index, 1);
      setSelectedLabels(updatedLabels);
    }
  };

  const getDatasetByLabel = (label) => {
    return [ ...generatedDatasets].find((dataset) => dataset.label === label);
  };

  const selectedDatasets = selectedLabels.map(getDatasetByLabel);

  return (
    <div className="HoxGeneratorApp">
      <header className="App-header-1">
        <h1>Hox-Generator</h1>
      </header>
      <DataGenerator onDataGenerated={handleDataGenerated} />
    </div>
  );
}

export default HoxGeneratorApp;
