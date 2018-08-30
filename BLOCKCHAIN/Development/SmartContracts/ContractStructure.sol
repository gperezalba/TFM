pragma solidity ^0.4.18;

//import "./mortal.sol";

contract ContractStructure {

    modifier onlyowner() {
        if (owner == msg.sender){
            _;
        } else {
            revert();
        }
    }

    uint someVar;
    address owner;

    mapping(address => Permission) myAddressMapping;
    
    struct Permission {
        bool isAllowed;
        uint maxAmount;
    }
    
    constructor() public payable {
        someVar = 0;
        owner = msg.sender;
        myAddressMapping[msg.sender] = Permission (true,5);
    }
    
    function setSomeVar(uint myVar) public onlyowner {
        someVar = myVar;
    }
    
    function getSomeVar() public view returns (uint) {
        return someVar;
    }
    
    function getContractBalance() public view returns(uint) {
        return address(this).balance;
    }

    function () public payable {
        
    }
}