import React, { useState, useRef, useEffect } from 'react';
import { ChromePicker } from 'react-color';

const LegendDropdown = ({ datasets, activeDatasets, toggleDataset, clearAll, handleColorChange, title }) => {
  const [searchTerm, setSearchTerm] = useState('');
  const [isOpen, setIsOpen] = useState(false);
  const [colorPickerLabel, setColorPickerLabel] = useState(null);
  const colorPickerRef = useRef(null);

  const filteredDatasets = datasets.filter(dataset =>
    dataset.label.toLowerCase().includes(searchTerm.toLowerCase())
  );

  const handleToggle = () => {
    setIsOpen(!isOpen);
  };

  const handleClearAll = () => {
    clearAll();
    setIsOpen(false); // Close the dropdown after clearing all
  };

  const handleColorPickerToggle = (label, event) => {
    event.stopPropagation(); // Prevent the dropdown from closing
    setColorPickerLabel(colorPickerLabel === label ? null : label);
  };

  // Function to set the selected color with alpha fixed to 0.5
  const setColorWithFixedAlpha = (color) => {
    return {
      ...color,
      rgba: {
        ...color.rgba,
        a: 0.5 // Set alpha to 0.5
      },
      hex: color.hex + '50' // Optional: append 80 for 50% opacity in hex
    };
  };

  // Close color picker if clicking outside
  useEffect(() => {
    const handleClickOutside = (event) => {
      if (colorPickerRef.current && !colorPickerRef.current.contains(event.target)) {
        setColorPickerLabel(null);
      }
    };
    document.addEventListener('mousedown', handleClickOutside);
    return () => {
      document.removeEventListener('mousedown', handleClickOutside);
    };
  }, [colorPickerRef]);

  return (
    <div className="legend-dropdown">
      <div className="dropdown-header" onClick={handleToggle}>
        <button className="dropdown-button">{title}</button>
        <button className="clear-all-button" onClick={handleClearAll}>Clear All</button>
      </div>
      <div className={`dropdown-content ${isOpen ? 'active' : ''}`}>
        <div className="search-container">
          <input
            type="text"
            placeholder="Search Motifs"
            value={searchTerm}
            onChange={e => setSearchTerm(e.target.value)}
            className="search-input"
            onClick={(e) => e.stopPropagation()} // Prevent dropdown from closing
          />
        </div>
        <div className="datasets-container">
          {filteredDatasets.map((dataset) => (
            <div key={dataset.label} className="dropdown-item">
              <div className="checkbox-container">
                <input
                  type="checkbox"
                  id={dataset.label}
                  checked={activeDatasets.includes(dataset.label)}
                  onChange={() => toggleDataset(dataset.label)}
                />
                <label htmlFor={dataset.label}>{dataset.label}</label>
              </div>
              <div className="color-picker-container" ref={colorPickerRef}>
                <button
                  className="color-picker-button"
                  style={{
                    backgroundColor: dataset.backgroundColor,
                    borderColor: dataset.backgroundColor, // Match border color to background color
                    borderWidth: '0px',
                    borderStyle: 'solid'
                  }}
                  onClick={(e) => handleColorPickerToggle(dataset.label, e)} // Prevent dropdown from closing
                >
                  <span className="color-picker-icon">▶</span>
                </button>
                {colorPickerLabel === dataset.label && (
                  <div className="color-picker-menu">
                    <ChromePicker
                      color={dataset.backgroundColor} // Allow full color selection
                      onChangeComplete={(color) => {
                        const updatedColor = setColorWithFixedAlpha(color);
                        handleColorChange(
                          datasets.findIndex(d => d.label === dataset.label),
                          updatedColor // Always update color with fixed alpha
                        );
                      }}
                      onChange={(color) => {
                        const updatedColor = setColorWithFixedAlpha(color);
                        handleColorChange(
                          datasets.findIndex(d => d.label === dataset.label),
                          updatedColor // Always update color with fixed alpha in real-time
                        );
                      }}
                    />
                  </div>
                )}
              </div>
            </div>
          ))}
        </div>
      </div>
    </div>
  );
};

export default LegendDropdown;
