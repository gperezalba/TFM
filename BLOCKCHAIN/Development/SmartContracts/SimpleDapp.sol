pragma solidity ^0.4.24;

contract SimpleDapp {
    
    function setSomeVar(uint myVar) public{
        someVar = myVar;
    }
    
    function getSomeVar() public view returns (uint) {
        return someVar;
    }
    
    function setSomeVarTimesFour(uint myVar) public {
        setSomeVar(4 * myVar);
    }
    uint someVar;
}