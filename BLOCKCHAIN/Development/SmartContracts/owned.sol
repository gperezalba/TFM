pragma solidity ^0.4.18;

contract owned {
    address owner;

    modifier onlyowner() {
        if (owner == msg.sender){
            _;
        } else {
            revert();
        }
    }

    constructor() public {
        owner = msg.sender;
    }
}