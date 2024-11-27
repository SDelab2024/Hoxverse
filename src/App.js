import React, { Suspense, lazy, useState, useEffect } from 'react';
import { BrowserRouter as Router, Route, Routes, Link } from 'react-router-dom';
import './App.css';
import DNASpinner from './components/DNASpinner';
import HoxGeneratorApp from './HoxGeneratorApp.js'; // Update the path accordingly

const ScatterPlot = lazy(() => import('./components/ScatterPlot'));
const ScatterPlot_Ranked = lazy(() => import('./components/Scatterplot_Ranked'));

function HomePage() {
  useEffect(() => {
    const handleScroll = () => {
      if (window.scrollY > 0) {
        document.body.classList.add('scrolled');
      } else {
        document.body.classList.remove('scrolled');
      }
    };

    window.addEventListener('scroll', handleScroll);
    
    // Cleanup the event listener
    return () => {
      window.removeEventListener('scroll', handleScroll);
    };
  }, []);

  return (
    <div className="HomePage">
      <header className="App-header">
        <div className='Welcome'>
          <TypingText text="Weelcome to Protein Engineering Lab" speed={50} />
        </div>
      </header>
      <div className='Options'>
        <h1>Presenting HOXverse</h1>
        <h2>Select the type of scatter plot you want to view:</h2>
        <nav className="navigation-links">
          <Link to="/unranked">Unranked Protein Scatter Plot</Link>
          <Link to="/ranked">Ranked Protein Scatter Plot</Link>
          <Link to="/hox-generator">Hox Searcher</Link>
        </nav>
        <br /><br /><br /><br />
        <h4 className="help">
          Note**
          <div className='note'>Incase there is no output result then that motif is not present in the following protein. <br /> Thank You for using HOXverse.</div>
        </h4>
        <br /><br /><br /><br /><br /><br />
        <div className="about"> 
          <h2 h2 className="about-title">Our Team</h2>
            <div className="team-container">
              <div className="team-member">
                {/* <img src="path/to/image.jpg" alt="Team Member Name" className="team-photo" /> */}
                  <h3 className="team-name">Dr. Soumya De</h3>
                    <p className="team-role">Principal Investigator</p>
                      <p className="team-bio">Associate Professor, Department of Bioscience and Biotechnology, Indian Institute of Technology Kharagpur</p>
                  </div>
              <div className="team-member">
                {/* <img src="path/to/image.jpg" alt="Team Member Name" className="team-photo" /> */}
                  <h3 className="team-name">Dr. Srijeeb Karmakar</h3>
                    <p className="team-role">Co-Principal Investigator</p>
                      <p className="team-bio">National Post Doctoral Fellow, Department of Bioscience and Biotechnology, Indian Institute of Technology Kharagpur</p>
                  </div> 
              <div className="team-member">
                {/* <img src="path/to/image.jpg" alt="Team Member Name" className="team-photo" /> */}
                  <h3 className="team-name">Jishnu Manglam </h3>
                    <p className="team-role">Developer</p>
                      <p className="team-bio">Undergrad Student, Department of Bioscience and Biotechnology, Indian Institute of Technology Kharagpur</p>
                  </div> 
              <div className="team-member">
                {/* <img src="path/to/image.jpg" alt="Team Member Name" className="team-photo" /> */}
                  <h3 className="team-name">Krishna Kant</h3>
                    <p className="team-role">Developer</p>
                      <p className="team-bio">Undergrad Student, Department of Bioscience and Biotechnology, Indian Institute of Technology Kharagpur</p>
                  </div>           
                </div>
              </div>
              <h4 className="help">
          Note**
          <div className='note'>If you are using Hoxverse data in your paper, please cite us.<br /> Thank You for using HOXverse.</div>
        </h4>
      </div>
    </div>
  );
}

function TypingText({ text, speed }) {
  const [displayText, setDisplayText] = useState('');

  useEffect(() => {
    let index = 0;
    const intervalId = setInterval(() => {
      if (index < text.length) {
        setDisplayText((prev) => prev + text.charAt(index));
        index++;
      } else {
        clearInterval(intervalId);
      }
    }, speed);

    return () => clearInterval(intervalId);
  }, [text, speed]);

  return (
    <h1 style={{ 
      fontSize: '50px',  // Increased font size
      fontWeight: 'bold', 
      textAlign: 'center', 
      margin: '0 auto',
      position: 'absolute',
      top: '50%',
      left: '50%',
      transform: 'translate(-50%, -50%)',
      zIndex: 1 // To ensure it appears above other elements
    }}>
      {displayText}
    </h1>
  );
}


function App() {
  return (
    <Router>
      <div className="App">
        <Routes>
          <Route path="/" element={<HomePage />} />
          <Route
            path="/unranked"
            element={
              <Suspense fallback={<DNASpinner />}>
                <ScatterPlot activeDataset="Dataset 1" />
              </Suspense>
            }
          />
          <Route
            path="/ranked"
            element={
              <Suspense fallback={<DNASpinner />}>
                <ScatterPlot_Ranked activeDataset="Dataset 1" />
              </Suspense>
            }
          />
          <Route
            path="/hox-generator"
            element={<HoxGeneratorApp />} // Added route
          />
        </Routes>
      </div>
    </Router>
  );
}

export default App;
