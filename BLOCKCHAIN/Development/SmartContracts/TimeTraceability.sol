pragma solidity ^0.4.18;

import "./mortal.sol";

contract TimeTraceability is mortal {

    address public owner;

    struct Data {
        bytes32[] hashCGGTTS;
        bool active;
        mapping(bytes32 => bool) valids;
    }

    mapping(address => Data) members;

    constructor(){
        owner = msg.sender;
    }

    function addMember(address newAddress) public onlyowner returns (bool) {
        require(!members[newAddress].active, "Member is already active");
        members[newAddress].active = true;
        return true;
    }

    function updateData(bytes32 data) public returns (uint length){
        require(members[msg.sender].active, "Not a member");
        members[msg.sender].valids[data] = true;
        return members[msg.sender].hashCGGTTS.push(data);
    }

    function getData(uint id, address ownerAddress) public view returns (bytes32) {
        require(members[ownerAddress].active, "Not a member");
        return members[ownerAddress].hashCGGTTS[id];
    }

    function getDataArray(address ownerAddress) public view returns (bytes32[]) {
        require(members[ownerAddress].active, "Not a member");
        return members[ownerAddress].hashCGGTTS;
    }

    function validateHash(bytes32 data, address ownerAddress) public view returns (bool) {
        require(members[ownerAddress].active, "Not a member");
        return members[ownerAddress].valids[data];
    }

    function() payable {
        revert();
    }
}

